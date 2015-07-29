import os
import logging
from random import shuffle
import numpy as np

from sklearn.externals import joblib

from django.db import models
from django.db.models import Q
from django.contrib.postgres.fields import ArrayField
from django.conf import settings
from django.utils import timezone
from django.core.validators import MaxValueValidator, MinValueValidator

from progressbar import ProgressBar, Percentage, Bar, ETA
from gensim.models import Doc2Vec, Phrases

from .constants import MODEL_STATES, FIELDS_FOR_MODEL, LSHF_STATES
from .utils import paper2tokens, TaggedDocumentsIterator, MyLSHForest
from .exceptions import InvalidState

from core.models import TimeStampedModel
from core.utils import pad_vector
from library.models import Paper, Journal
from core.constants import TIME_LAPSE_CHOICES

logger = logging.getLogger(__name__)


class ModelManager(models.Manager):
    def create(self, **kwargs):
        model = super(ModelManager, self).create(**kwargs)

        # paper fields used to generate model data
        if 'fields' in kwargs:
            if not isinstance(kwargs['fields'], list):
                raise TypeError('<fields> must be list of field strings')
        fields = kwargs.get('fields', ['title', 'abstract'])
        for field in fields:
            FieldUseInModel.objects.create(model=model, field=field)

        # Init doc2vec from gensim
        model.init_doc2vec()

        return model

    def get(self, **kwargs):
        print('Warning: Not loading Doc2Vec. Use load() instead')
        return super(ModelManager, self).get(**kwargs)

    def load(self, **kwargs):
        model = super(ModelManager, self).get(**kwargs)
        try:
            model.doc2vec = Doc2Vec.load(os.path.join(settings.NLP_DOC2VEC_PATH,
                                                      '{0}.mod'.format(
                                                          model.name)))
        except FileNotFoundError:
            raise FileNotFoundError
        return model


class Model(TimeStampedModel):
    """
    Model for a doc2vec type of Natural Language Modeling leveraging gensim

    Example of use case for training a new model:
    1) Create model
    >>> model = Model.objects.create(name='test')
    2) Prepare and dump data.
    >>> papers = Paper.objects.all()
    >>> model.dump(papers)
    3) Build phraser and vocab (phraser: join common n-gram word. ie new_york,
    optional)
    >>> docs = model.load_documents()
    >>> phraser = model.build_phraser(docs, min_count=2)
    >>> docs = model.load_documents(phraser=phraser)
    >>> model.build_vocab(docs)
    4) Train, save and set_active
    >>> model.train(docs, passes=10, shuffle_=True)
    >>> model.save()
    >>> model.set_active()
    5) Populate library is needed
    >>> model.save_journal_vec_from_bulk()
    >>> model.save_paper_vec_from_bulk()
    6) Start embedding new paper <paper> as they come in
    >>> model.infer_paper(paper)
    """

    name = models.CharField(max_length=128, blank=False, null=False,
                            unique=True)

    # when trained, model becomes active
    is_active = models.BooleanField(default=False)

    # doc2vec instance from gensim
    doc2vec = Doc2Vec()

    objects = ModelManager()

    # Model parameters
    # ----------------------
    # `dm` defines the training algorithm. By default (`dm=1`), 'distributed
    # memory' (PV-DM) is used. Otherwise, `distributed bag of words` (PV-DBOW)
    # is employed.
    dm = models.IntegerField(default=1)

    # `size` is the dimensionality of the feature vectors.
    size = models.IntegerField(default=100,
           validators=[MaxValueValidator(settings.NLP_MAX_VECTOR_SIZE)])

    # `window` is the maximum distance between the predicted word and context
    # words used for prediction within a document.
    window = models.IntegerField(default=8)

    # `alpha` is the initial learning rate (will linearly drop to zero as
    # training progresses).
    alpha = models.FloatField(default=0.025)

    # `seed` = for the random number generator. Only runs with a single worker
    # will be deterministically reproducible because of the ordering randomness
    # in multi-threaded runs.
    seed = models.IntegerField(default=0)

    # `min_count` = ignore all words with total frequency lower than this.
    min_count = models.IntegerField(default=2)

    # `max_vocab_size` = limit RAM during vocabulary building; if there are more
    #  unique words than this, then prune the infrequent ones. Every 10 million
    # word types need about 1GB of RAM. Set to `None` for no limit (default).
    max_vocab_size = models.IntegerField(default=None, null=True, blank=True)

    # `sample` = threshold for configuring which higher-frequency words are
    # randomly downsampled; default is 0 (off), useful value is 1e-5.
    sample = models.FloatField(default=0.)

    # `workers` = use this many worker threads to train the model
    # (=faster training with multicore machines).
    workers = models.IntegerField(default=1)

    # `hs` = if 1 (default), hierarchical sampling will be used for model
    # training (else set to 0).
    hs = models.IntegerField(default=1)

    # `negative` = if > 0, negative sampling will be used, the int for negative
    # specifies how many "noise words" should be drawn (usually between 5-20).
    negative = models.IntegerField(default=0)

    # `dm_mean` = if 0 (default), use the sum of the context word vectors.
    # If 1, use the mean. Only applies when dm is used in non-concatenative mode.
    dm_mean = models.IntegerField(default=0)

    # `dm_concat` = if 1, use concatenation of context vectors rather than
    # sum/average; default is 0 (off). Note concatenation results in a
    # much-larger model, as the input is no longer the size of one (sampled or
    # arithmatically combined) word vector, but the size of the tag(s) and all
    # words in the context strung together.
    dm_concat = models.IntegerField(default=0)

    # `dm_tag_count` = expected constant number of document tags per document,
    # when using dm_concat mode; default is 1.
    dm_tag_count = models.IntegerField(default=1)

    # `dbow_words` if set to 1 trains word-vectors (in skip-gram fashion)
    # simultaneous with DBOW doc-vector training; default is 0 (faster training
    # of doc-vectors only).
    dbow_words = models.IntegerField(default=0)

    # use in model clean
    model_arguments = [
        ('dm', 'dm'),
        ('size', 'vector_size'),
        ('window', 'window'),
        ('alpha', 'alpha'),
        ('seed', 'seed'),
        ('min_count', 'min_count'),
        ('max_vocab_size', 'max_vocab_size'),
        ('sample', 'sample'),
        ('workers', 'workers'),
        ('hs', 'hs'),
        ('negative', 'negative'),
        ('dm_mean', 'cbow_mean'),
        ('dm_concat', 'dm_concat'),
        ('dm_tag_count', 'dm_tag_count'),
        ('dbow_words', 'dbow_words'),
    ]

    def __str__(self):
        return self.name

    def init_doc2vec(self):
        """Instantiate new doc2vec with model field attributes
        """
        self.doc2vec = Doc2Vec(
            dm=self.dm,
            size=self.size,
            window=self.window,
            alpha=self.alpha,
            seed=self.seed,
            min_count=self.min_count,
            max_vocab_size=self.max_vocab_size,
            sample=self.sample,
            workers=self.workers,
            hs=self.hs,
            negative=self.negative,
            dm_mean=self.dm_mean,
            dm_concat=self.dm_concat,
            dm_tag_count=self.dm_tag_count,
            dbow_words=self.dbow_words
        )

    def save_db_only(self, *args, **kwargs):
        self.full_clean()
        super(Model, self).save(*args, **kwargs)

    def save(self, *args, **kwargs):
        self.full_clean()
        super(Model, self).save(*args, **kwargs)
        if not os.path.exists(settings.NLP_DOC2VEC_PATH):
            os.makedirs(settings.NLP_DOC2VEC_PATH)
        self.doc2vec.save(
            os.path.join(settings.NLP_DOC2VEC_PATH, '{0}.mod'.format(self.name)))

    def delete(self, *args, **kwargs):
        try:
            os.remove(
                os.path.join(settings.NLP_DOC2VEC_PATH, '{0}.mod'.format(self.name)))
        except FileNotFoundError:
            pass
        super(Model, self).delete(*args, **kwargs)

    def clean(self):
        """Check  model field attributes and doc2vec attributes are in agreements
        """
        for attr in self.model_arguments:
            if not getattr(self, attr[0]) == getattr(self.doc2vec, attr[1]):
                raise ValueError('{attr} does not match'.format(attr=attr[0]))

    def dump(self, papers):
        """Dump papers data to pre-process text files.
        Papers data are separated in multiple file to spare memory when building
        model. File are composed of N (CHUNK_SIZE) documents. Each line of the
        files start with the primary key of the document and then the pre-process
        string corresponding to the document field(s) dumped.
         ex:
            123: the title of the document #1 . the abstract of the document .
            124: the title of the document #2 . the abstract of the document .

        Dump paper date to files.
        Note: It is useful to order paper randomly (order_by('?')) to avoid bias
        """

        fields = list(self.paperfields.all().values_list('field', flat='True'))
        tot = papers.count()
        file_count = 0
        file = None

        if not os.path.exists(settings.NLP_DATA_PATH):
            os.makedirs(settings.NLP_DATA_PATH)

        sub_update_step = 20
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=int(tot/sub_update_step),
                           redirect_stderr=True).start()

        for count, paper in enumerate(papers):

            if not count % settings.NLP_CHUNK_SIZE:
                if file:
                    file.close()
                file = open(os.path.join(settings.NLP_DATA_PATH,
                                         '{0:06d}.txt'.format(file_count)),
                            'w+')
                file_count += 1

            # write header
            if paper.journal:
                j_pk = paper.journal.pk
            else:
                j_pk = 0
            file.write('{pk}, j_{j_pk}: '.format(pk=paper.pk, j_pk=j_pk))

            # line body
            line_val = paper2tokens(paper, fields=fields)

            # write to file
            file.write(' '.join(line_val).strip())

            # write new line
            file.write('\n')

            # update progress bar
            if not count % sub_update_step:
                pbar.update(count/sub_update_step)
        # close progress bar
        pbar.finish()
        file.close()

    def load_documents(self, phraser=None):
        """Load text files stored in data_path to documents list
        """
        return TaggedDocumentsIterator(settings.NLP_DATA_PATH, phraser=phraser)

    @staticmethod
    def build_phraser(documents, min_count=2, threshold=10.0):
        """build phraser
        """
        phraser = Phrases(map(lambda doc: doc.words, documents),
                          min_count=min_count, threshold=threshold)
        return phraser

    def build_vocab(self, documents):
        """build vocabulary
        """
        self.doc2vec.build_vocab(documents)

    def train(self, documents, passes=10, shuffle_=True, alpha=0.025,
              min_alpha=0.001):
        """Train model

        Using explicit multiple-pass, alpha-reduction approach as sketched in
        gensim doc2vec blog post (http://rare-technologies.com/doc2vec-tutorial)
        """

        # Init
        alpha_delta = (alpha - min_alpha) / passes
        if shuffle_:  # load all doc to memory
            documents = [doc for doc in documents]

        # train
        for epoch in range(passes):
            print('Epoch {epoch}/{passes}\n'.format(epoch=epoch, passes=passes))
            if shuffle_:
                # shuffling before training is in general a good habit
                shuffle(documents)
            self.alpha = alpha
            self.doc2vec.alpha = alpha
            self.doc2vec.min_alpha = alpha
            self.doc2vec.train(documents)
            alpha -= alpha_delta

    def save_journal_vec_from_bulk(self):
        """Store inferred journal vector from training in db
        """
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(self.doc2vec.docvecs.doctags),
                           redirect_stderr=True).start()
        count = 0
        # TODO: replace with bulk update / create ?
        for pk in self.doc2vec.docvecs.doctags:
            if pk.startwith('j') and not pk == 'j_0':
                try:
                    journal = Journal.objects.get(pk=int(pk[2:]))
                    vector = self.doc2vec.docvecs[pk].tolist()
                    pv, _ = JournalVectors.objects.get_or_create(model=self,
                                                                 journal=journal)
                    pv.vector = vector
                    pv.save()
                except Journal.DoesNotExist:
                    print('Journal {pk} did not match'.format(pk=pk))
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

    def save_paper_vec_from_bulk(self):
        """Store inferred paper vector from training in db
        """
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(self.doc2vec.docvecs.doctags),
                           redirect_stderr=True).start()
        count = 0
        # TODO: replace with bulk update / create ?
        for pk in self.doc2vec.docvecs.doctags:
            if not pk.startwith('j'):
                try:
                    paper = Paper.objects.get(pk=int(pk))
                    vector = self.doc2vec.docvecs[pk].tolist()
                    pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                               paper=paper)
                    pv.set_vector(vector)
                    pv.save()
                except Paper.DoesNotExist:
                    print('Paper {pk} did not match'.format(pk=pk))
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

    def save_all_vec_from_bulk(self):
        """Store inferred paper and journal vectors from training in db
        """
        self.save_paper_vec_from_bulk()
        self.save_journal_vec_from_bulk()

    def infer_paper(self, paper_pk, alpha=0.05, min_alpha=0.001, passes=5):
        """Infer model vector for paper
        """
        fields = list(self.paperfields.all().values_list('field', flat='True'))

        pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                   paper_id=paper_pk)
        paper = Paper.objects.get(id=paper_pk)
        doc_words = paper2tokens(paper, fields=fields)

        # NB: infer_vector user explicit alpha-reduction with multi-passes
        vector = self.doc2vec.infer_vector(doc_words,
                                           alpha=alpha,
                                           min_alpha=min_alpha,
                                           steps=passes)
        pv.set_vector(vector.tolist())
        pv.save()

        return vector

    def infer_papers(self, papers, **kwargs):
        """Infer model vector for papers
        """

        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=papers.count(), redirect_stderr=True).start()
        count = 0
        for paper in papers:
            self.infer_paper(paper, **kwargs)
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

    def set_active(self):
        self.is_active = True
        self.save_db_only()

    class Meta:
        ordering = ['name', ]

    @classmethod
    def infer_paper_all_models(cls, paper, **kwargs):
        """Infer all model vector for paper
        """
        model_names = list(cls.objects.all().values_list('name', flat='True'))
        for model_name in model_names:
            model = cls.objects.get(name=model_name)
            model.infer_paper(paper, **kwargs)

    @classmethod
    def infer_papers_all_models(cls, papers, **kwargs):
        """Infer all model vector for all papers
        """
        model_names = list(cls.objects.all().values_list('name', flat='True'))
        for model_name in model_names:
            model = cls.objects.get(name=model_name)
            model.infer_papers(papers, **kwargs)


class FieldUseInModel(TimeStampedModel):
    """Store which fields of paper and related data are used in model training
    """
    model = models.ForeignKey(Model, related_name='paperfields')

    field = models.CharField(max_length=100, choices=FIELDS_FOR_MODEL)

    class Meta:
        unique_together = [('model', 'field')]

    def __str__(self):
        return '{field}'.format(field=self.field)


class PaperVectorsManager(models.Manager):

    def create(self, **kwargs):
        # enforce that model is active
        if not kwargs.get('model').is_active:
            raise ValueError('model {0} is not active'
                             .format(kwargs.get('model').name))

        return super(PaperVectorsManager, self).create(**kwargs)


class PaperVectors(TimeStampedModel):
    """Paper - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None

    Use set_vector() to pad and set vector list
    """

    paper = models.ForeignKey(Paper, related_name='vectors')

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    is_in_full_lsh = models.BooleanField(default=False)

    objects = PaperVectorsManager()

    class meta:
        unique_together = ('paper', 'model')

    def __str__(self):
        return '{short_title}/{name}'.format(short_title=self.paper.short_title,
                                             name=self.model.name)

    def infer_paper(self, **kwargs):
        vector = self.model.infer_paper(self.paper, **kwargs)
        self.set_vector(vector)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()


class PaperNeighbors(TimeStampedModel):

    model = models.ForeignKey(Model)

    paper = models.ForeignKey(Paper)

    neighbor1 = models.ForeignKey(Paper, related_name='neighbor1', null=True)

    neighbor2 = models.ForeignKey(Paper, related_name='neighbor2', null=True)

    neighbor3 = models.ForeignKey(Paper, related_name='neighbor3', null=True)

    neighbor4 = models.ForeignKey(Paper, related_name='neighbor4', null=True)

    neighbor5 = models.ForeignKey(Paper, related_name='neighbor5', null=True)


class JournalVectors(TimeStampedModel):
    """Journal - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None.

    Use set_vector() to pad and set vector list
    """

    journal = models.ForeignKey(Journal)

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class meta:
        unique_together = ('journal', 'model')

    def __str__(self):
        return '{short_title}/{name}'\
            .format(short_title=self.journal.short_title, name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()


class LSHManager(models.Manager):

    def create(self, **kwargs):
        model = super(LSHManager, self).create(**kwargs)
        # Init lsh
        model.full_update(**kwargs)
        return model

    def get(self, **kwargs):
        # print('Warning: Not loading LSH object. Use load() instead')
        return super(LSHManager, self).get(**kwargs)

    def load(self, **kwargs):
        obj = super(LSHManager, self).get(**kwargs)
        try:
            obj.lsh = joblib.load(os.path.join(settings.NLP_LSH_PATH,
                                       '{0}.lsh'.format(obj.model.name)))
        except FileNotFoundError:
            raise FileNotFoundError
        return obj


class LSH(TimeStampedModel):
    """Local Sensitive Hashing to retrieve approximate k-neighbors
    """

    model = models.OneToOneField(Model, related_name='lsh')

    state = models.CharField(default='', choices=LSHF_STATES, max_length=3)

    time_lapse = models.IntegerField(default=None,
                                     null=True, blank=True,
                                     choices=TIME_LAPSE_CHOICES,
                                     verbose_name='In the past for')

    lsh = MyLSHForest()

    objects = LSHManager()

    def save(self, *args, **kwargs):
        if not os.path.exists(settings.NLP_LSH_PATH):
            os.makedirs(settings.NLP_LSHF_PATH)
        joblib.dump(self.lsh, os.path.join(settings.NLP_LSH_PATH,
                                       '{0}.lsh'.format(self.model.name)))
        # Model is now Idle
        self.state = 'IDL'
        super(LSH, self).save(*args, **kwargs)

    def save_db_only(self, *args, **kwargs):
        super(LSH, self).save(*args, **kwargs)

    def full_update(self, n_estimators=10, n_candidates=50):

        # Register status as busy
        self.state = 'BUS'
        self.save_db_only()

        # Init LSHF
        self.lsh = MyLSHForest(n_estimators=n_estimators,
                               n_candidates=n_candidates)

        # Get data
        if self.time_lapse:   # time range
            from_date = (timezone.now() -
                         timezone.timedelta(
                             days=self.time_lapse)).date()

            data = PaperVectors.objects\
                            .filter(Q(model=self.model) &
                                    (Q(paper__date_ep__gt=from_date) |
                                     (Q(paper__date_pp__gt=from_date) &
                                      Q(paper__date_ep=None))))\
                            .values_list('pk', 'paper__pk', 'vector')
        else:
            data = PaperVectors.objects\
                .filter(model=self.model)\
                .values_list('pk', 'paper__pk', 'vector')
        vec_size = self.model.size

        # Reshape data
        pv_pks = []
        X = np.zeros((data.count(), vec_size))
        for i, dat in enumerate(data):
            # store PaperVector pk
            pv_pks.append(data[i][0])
            # store paper pk
            self.lsh.pks.append(data[i][1])
            # build input matrix for fit
            X[i, :] = data[i][2][:vec_size]

        # Train
        self.lsh.fit(X)

        # Save
        self.save()

        # Update PaperVector.is_in_full_lsh
        if not self.time_lapse:  # this is a the full LSH
            PaperVectors.objects.filter(pk__in=pv_pks).update(is_in_full_lsh=True)

    def partial_update(self):
        """Only partial update LSH with new incoming papers (i.e. does not remove
        old ones, valid only if self.time_lapse=None)

        NB: we can not do partial on time_lapse LSH because paper cannot be removed from
        training set, only added.
        """

        if not self.time_lapse:
            # Register status as busy
            self.state = 'BUS'
            self.save_db_only()

            # Get data
            data = PaperVectors.objects\
                .filter(model=self.model, is_in_full_lsh=False)\
                .values_list('pk', 'paper__pk', 'vector')
            vec_size = self.model.size

            # Reshape data
            pv_pks = []
            X = np.zeros((data.count(), vec_size))
            for i, dat in enumerate(data):
                # store PaperVector pk
                pv_pks.append(data[i][0])
                # store paper pk
                self.lsh.pks.append(data[i][1])
                # build input matrix for fit
                X[i, :] = data[i][2][:vec_size]

            # Train
            self.lsh.partial_fit(X)

            # Save
            self.save()

            # Update PaperVector.is_in_full_lsh
            PaperVectors.objects.filter(pk__in=pv_pks)\
                .update(is_in_full_lsh=True)
        else:
            logger.warning('LSH ({id}) attribute time_lapse is not null, '
                           'calling full_update() instead'.format(id=self.id))
            self.full_update()

    def kneighbors(self, vec, n_neighbors=10000):

        if self.state == 'IDL':
            if isinstance(vec, list):
                vec = np.array(vec)
            distances, indices = self.lsh.kneighbors(vec[:self.model.size],
                                                     n_neighbors=n_neighbors,
                                                     return_distance=True)
            # convert indices to paper pk
            for i, index in enumerate(indices):
                for j, k in enumerate(index):
                    indices[i][j] = self.lsh.pks[k]
            return indices
        else:
            raise InvalidState('LSH is not Idle. current state is {0}'
                               .format(self.state))

    def mono_task(self, *args, **kwargs):
        """Use form tasks.py
        """

        try:
            task = kwargs['task']
        except KeyError:
            raise KeyError('key task must defined')

        if task == 'update':
            self.full_update()
            return 0
        elif task == 'partial_update':
            self.partial_update()
            return 0
        elif task == 'kneighbors':
            try:
                vec = kwargs['vec']
                n_neighbors = kwargs['n_neighbors']
            except ValueError as e:
                raise ValueError(e)
            pks = self.kneighbors(vec, n_neighbors)
            return pks
        else:
            print(task)
            raise ValueError('Unknown task action')