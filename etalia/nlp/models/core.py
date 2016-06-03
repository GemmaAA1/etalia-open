# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import logging
import glob
from random import shuffle
import collections

import numpy as np
from gensim import matutils
from sklearn.externals import joblib
from progressbar import ProgressBar, Percentage, Bar, ETA
from django.db import models
from django.conf import settings
from django.db.models import QuerySet
from django.utils import timezone
from django.core.validators import MaxValueValidator
from django.core.exceptions import ValidationError
from gensim.models import Doc2Vec

from gensim.models import Phrases

from etalia.core.models import TimeStampedModel

from etalia.library.models import Paper, Journal, AuthorPaper
from etalia.threads.models import Thread, ThreadPost, ThreadComment

from .library import ModelLibraryMixin, PaperVectors, PaperNeighbors, \
    JournalVectors, PaperEngineScoringMixin
from .threads import ModelThreadMixin, ThreadNeighbors, ThreadVectors

from ..constants import FIELDS_FOR_MODEL
from ..utils import obj2tokens, TaggedDocumentsIterator, model_attr_getter, \
    model_attr_setter
from ..mixins import S3Mixin

logger = logging.getLogger(__name__)


class ModelManager(models.Manager):
    def create(self, **kwargs):
        # starting popping text_fields key if any
        default_text_fields = [field[0] for field in FIELDS_FOR_MODEL]
        text_fields = kwargs.pop('text_fields', default_text_fields)

        # instantiate object and save to db
        obj = self.model(**kwargs)
        obj.save_db_only(using=self.db)

        # Add text_fields used to generate model data
        # First test validity
        if not isinstance(text_fields, list):
            raise ValidationError('<text_fields> must be list of field strings')
        if not set(text_fields).issubset(default_text_fields):
            raise ValidationError(
                '<text_fields> not subset of FIELDS_FOR_MODEL')

        for text_field in text_fields:
            fuim, _ = TextField.objects.get_or_create(text_field=text_field)
            obj.text_fields.add(fuim)

        obj.full_clean()
        obj.save_db_only(using=self.db)

        # Init doc2vec from gensim
        obj.init_doc2vec()

        obj.save(using=self.db)

        return obj

    def load(self, **kwargs):
        """Get model instance and load model data"""
        obj = super(ModelManager, self).get(**kwargs)
        try:
            # if not on volume try download from s3
            if not os.path.isfile(os.path.join(settings.NLP_MODELS_PATH,
                                               '{}.mod'.format(obj.name))):
                if getattr(settings, 'NLP_MODELS_BUCKET_NAME', ''):
                    obj.pull_from_s3()
            obj._doc2vec = Doc2Vec.load(os.path.join(
                settings.NLP_MODELS_PATH, '{0}.mod'.format(obj.name)))
        except EnvironmentError:  # OSError or IOError...
            raise
        return obj


class Model(ModelThreadMixin,
            ModelLibraryMixin,
            S3Mixin,
            TimeStampedModel):
    """Natural Language Processing Class based on Doc2Vec from Gensim
    """
    # For S3 Mixin
    BUCKET_NAME = getattr(settings, 'NLP_MODELS_BUCKET_NAME', '')
    PATH = getattr(settings, 'NLP_MODELS_PATH', '')
    AWS_ACCESS_KEY_ID = getattr(settings, 'AWS_ACCESS_KEY_ID', '')
    AWS_SECRET_ACCESS_KEY = getattr(settings, 'AWS_SECRET_ACCESS_KEY', '')

    name = models.CharField(max_length=128, blank=False, null=False,
                            unique=True)

    # when trained, model becomes active. If False model is not load as a model
    # in the workers
    is_active = models.BooleanField(default=False)

    # doc2vec instance from gensim
    _doc2vec = Doc2Vec()

    # fields of Paper used as text inputs during model training
    text_fields = models.ManyToManyField('TextField',
                                         related_name='text_fields')

    upload_state = models.CharField(max_length=3,
                                    choices=(('IDL', 'Idle'),
                                             ('ING', 'Uploading')),
                                    default='IDL')

    objects = ModelManager()

    # Model parameters (mirroring Doc2Vec attributes)
    # ----------------------
    # `dm` defines the training algorithm. By default (`dm=1`), 'distributed
    # memory' (PV-DM) is used. Otherwise, `distributed bag of words` (PV-DBOW)
    # is employed.
    _dm = models.IntegerField(default=1)

    # `size` is the dimensionality of the feature vectors.
    _size = models.IntegerField(default=128,
                                validators=[MaxValueValidator(
                                    settings.NLP_MAX_VECTOR_SIZE)])

    # `window` is the maximum distance between the predicted word and context
    # words used for prediction within a document.
    _window = models.IntegerField(default=8)

    # `alpha` is the initial learning rate (will linearly drop to zero as
    # training progresses).
    _alpha = models.FloatField(default=0.025)

    # `seed` = for the random number generator. Only runs with a single worker
    # will be deterministically reproducible because of the ordering randomness
    # in multi-threaded runs.
    _seed = models.IntegerField(default=0)

    # `min_count` = ignore all words with total frequency lower than this.
    _min_count = models.IntegerField(default=2)

    # `max_vocab_size` = limit RAM during vocabulary building; if there are more
    #  unique words than this, then prune the infrequent ones. Every 10 million
    # word types need about 1GB of RAM. Set to `None` for no limit (default).
    _max_vocab_size = models.IntegerField(default=None, null=True, blank=True)

    # `sample` = threshold for configuring which higher-frequency words are
    # randomly downsampled; default is 0 (off), useful value is 1e-5.
    _sample = models.FloatField(default=0.)

    # `workers` = use this many worker threads to train the model
    # (=faster training with multicore machines).
    _workers = models.IntegerField(default=1)

    # `hs` = if 1 (default), hierarchical sampling will be used for model
    # training (else set to 0).
    _hs = models.IntegerField(default=1)

    # `negative` = if > 0, negative sampling will be used, the int for negative
    # specifies how many "noise words" should be drawn (usually between 5-20).
    _negative = models.IntegerField(default=0)

    # `dm_mean` = if 0 (default), use the sum of the context word vectors.
    # If 1, use the mean. Only applies when dm is used in non-concatenative mode.
    _dm_mean = models.IntegerField(default=0)

    # `dm_concat` = if 1, use concatenation of context vectors rather than
    # sum/average; default is 0 (off). Note concatenation results in a
    # much-larger model, as the input is no longer the size of one (sampled or
    # arithmatically combined) word vector, but the size of the tag(s) and all
    # words in the context strung together.
    _dm_concat = models.IntegerField(default=0)

    # `dm_tag_count` = expected constant number of document tags per document,
    # when using dm_concat mode; default is 1.
    _dm_tag_count = models.IntegerField(default=1)

    # `dbow_words` if set to 1 trains word-vectors (in skip-gram fashion)
    # simultaneous with DBOW doc-vector training; default is 0 (faster training
    # of doc-vectors only).
    _dbow_words = models.IntegerField(default=0)

    # corresponding setter/getter to sync model attribute and Doc2Vec attributes
    dm = property(fset=model_attr_setter('_dm', 'dm'),
                  fget=model_attr_getter('_dm'))
    size = property(fset=model_attr_setter('_size', 'vector_size'),
                    fget=model_attr_getter('_size'))
    window = property(fset=model_attr_setter('_window', 'window'),
                      fget=model_attr_getter('_window'))
    alpha = property(fset=model_attr_setter('_alpha', 'alpha'),
                     fget=model_attr_getter('_alpha'))
    seed = property(fset=model_attr_setter('_seed', 'seed'),
                    fget=model_attr_getter('_seed'))
    min_count = property(fset=model_attr_setter('_min_count', 'min_count'),
                         fget=model_attr_getter('_min_count'))
    max_vocab_size = property(fset=model_attr_setter('_max_vocab_size',
                                                     'max_vocab_size'),
                              fget=model_attr_getter('_max_vocab_size'))
    sample = property(fset=model_attr_setter('_sample', 'sample'),
                      fget=model_attr_getter('_sample'))
    workers = property(fset=model_attr_setter('_workers', 'workers'),
                       fget=model_attr_getter('_workers'))
    hs = property(fset=model_attr_setter('_hs', 'hs'),
                  fget=model_attr_getter('_hs'))
    negative = property(fset=model_attr_setter('_negative', 'negative'),
                        fget=model_attr_getter('_negative'))
    dm_mean = property(fset=model_attr_setter('_dm_mean', 'cbow_mean'),
                       fget=model_attr_getter('_dm_mean'))
    dm_concat = property(fset=model_attr_setter('_dm_concat', 'dm_concat'),
                         fget=model_attr_getter('_dm_concat'))
    dm_tag_count = property(fset=model_attr_setter('_dm_tag_count',
                                                   'dm_tag_count'),
                            fget=model_attr_getter('_dm_tag_count'))
    dbow_words = property(fset=model_attr_setter('_dbow_words', 'dbow_words'),
                          fget=model_attr_getter('_dbow_words'))

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name', ]

    def init_doc2vec(self):
        """Instantiate new doc2vec with model attributes
        """
        self._doc2vec = Doc2Vec(
            dm=self._dm,
            size=self._size,
            window=self._window,
            alpha=self._alpha,
            seed=self._seed,
            min_count=self._min_count,
            max_vocab_size=self._max_vocab_size,
            sample=self._sample,
            workers=self._workers,
            hs=self._hs,
            negative=self._negative,
            dm_mean=self._dm_mean,
            dm_concat=self._dm_concat,
            dm_tag_count=self._dm_tag_count,
            dbow_words=self._dbow_words
        )

    def save_db_only(self, *args, **kwargs):
        super(Model, self).save(*args, **kwargs)

    def save(self, *args, **kwargs):
        # save files to local volume
        if not os.path.exists(settings.NLP_MODELS_PATH):
            os.makedirs(settings.NLP_MODELS_PATH)
        self._doc2vec.save(
            os.path.join(settings.NLP_MODELS_PATH,
                         '{0}.mod'.format(self.name)))
        # push files to s3
        if self.BUCKET_NAME:
            self.upload_state = 'ING'
            self.save_db_only()
            self.push_to_s3(ext='mod')
            self.upload_state = 'IDL'
        # save to db
        self.save_db_only(*args, **kwargs)

    def delete(self, s3=False, **kwargs):
        try:
            for filename in glob.glob(os.path.join(settings.NLP_MODELS_PATH,
                                                   '{0}.mod*'.format(
                                                       self.name))):
                os.remove(filename)
        except IOError:
            pass

        # delete on amazon s3
        if s3:
            try:
                if self.BUCKET_NAME:
                    self.delete_on_s3()
            except IOError:
                pass
        super(Model, self).delete(**kwargs)

    def dump(self, papers, data_path=None):
        """Dump matches data to pre-process text files.
        Papers data are separated in multiple file to spare memory when building
        model. File are composed of N (CHUNK_SIZE) documents. Each line of the
        files start with the primary key of the document and then the
        pre-process string corresponding to the document field(s) dumped.
         ex:
            123: the title of the document #1 . the abstract of the document .
            124: the title of the document #2 . the abstract of the document .

        Dump paper date to files.
        Note: It is useful to order paper randomly (order_by('?')) to avoid bias
        """

        text_fields = [tx.text_field for tx in self.text_fields.all()]
        if isinstance(papers, QuerySet):
            tot = papers.count()
        else:
            tot = len(papers)
        file_count = 0
        fid = None

        if not data_path:
            data_path = settings.NLP_DATA_PATH
        if not os.path.exists(data_path):
            os.makedirs(data_path)

        sub_update_step = 20
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=int(tot / sub_update_step),
                           redirect_stderr=True).start()

        for count, paper in enumerate(papers):

            if not paper.abstract or not paper.is_trusted:
                logger.warning('DISCARDING paper {pk} (is_trusted=False or '
                               'empty abstract)'.format(pk=paper.id))
                continue

            # rotate file when # lines is NLP_CHUNK_SIZE
            if not count % settings.NLP_CHUNK_SIZE or not fid:
                if fid:
                    fid.close()
                fid = open(os.path.join(data_path,
                                        '{0:06d}.txt'.format(file_count)), 'w+')
                file_count += 1

            # write header
            if paper.journal:
                j_pk = paper.journal.pk
            else:
                j_pk = 0
            # fid.write('{pk}, j_{j_pk}: '.format(pk=paper.pk, j_pk=j_pk).encode('utf-8'))
            fid.write('{pk}, j_{j_pk}: '.format(pk=paper.pk, j_pk=j_pk))

            # line body
            line_val = obj2tokens(paper, fields=text_fields)

            # write to file
            # fid.write(u' '.join(line_val).strip().encode('utf-8'))
            fid.write(u' '.join(line_val).strip())

            # write new line
            fid.write('\n')

            # update progress bar
            if not count % sub_update_step:
                pbar.update(count / sub_update_step)
        # close progress bar
        pbar.finish()
        if fid:
            fid.close()

    @staticmethod
    def load_documents(data_path=None, phraser=None, ):
        """Load text files stored in data_path to documents list
        """
        if not data_path:
            data_path = settings.NLP_DATA_PATH
        if phraser:
            return TaggedDocumentsIterator(data_path, phraser=phraser)
        else:
            return TaggedDocumentsIterator(data_path)

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
        self.deactivate()

        self._doc2vec.build_vocab(documents)

    def train(self, documents, passes=10, shuffle_=True, alpha=0.025,
              min_alpha=0.001):
        """Train model

        Using explicit multiple-pass, alpha-reduction approach as sketched in
        gensim doc2vec blog post (http://rare-technologies.com/doc2vec-tutorial)
        """
        self.deactivate()

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
            self._doc2vec.alpha = alpha
            self._doc2vec.min_alpha = alpha
            self._doc2vec.train(documents)
            alpha -= alpha_delta
        self.save()

    def build_vocab_and_train(self, **kwargs):
        """Build vocabulary (with phraser is set) and train model"""
        # Build vocabulary and phraser
        docs = self.load_documents()
        phraser = self.build_phraser(docs, **kwargs)
        docs = self.load_documents(phraser=phraser)
        self.build_vocab(docs)
        # Train, save
        self.train(docs, **kwargs)

    def propagate(self):
        """Save Journal and Paper vectors and init MostSimilar
        (To be run after model training)"""
        logging.info('{} Populate journals...'.format(self.name))
        self.save_journal_vec_from_bulk()
        logging.info('{} Populate matches...'.format(self.name))
        self.save_paper_vec_from_bulk()
        logging.info('{} Build MostSimilar...'.format(self.name))

    def save_journal_vec_from_bulk(self):
        """Store inferred journal vector from training in db
        """
        self.check_active()

        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(self._doc2vec.docvecs.doctags),
                           redirect_stderr=True).start()
        count = 0
        for pk in self._doc2vec.docvecs.doctags:
            if pk.startswith('j') and not pk == 'j_0':
                try:
                    journal = Journal.objects.get(pk=int(pk[2:]))
                    vector = self._doc2vec.docvecs[pk]
                    # normalize
                    vector /= np.linalg.norm(vector)
                    # store
                    pv, _ = JournalVectors.objects.get_or_create(model=self,
                                                                 journal=journal)
                    pv.set_vector(vector)
                except Journal.DoesNotExist:
                    print('Journal {pk} did not match'.format(pk=pk))
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

    def save_paper_vec_from_bulk(self):
        """Store inferred paper vector from training in db
        """
        self.check_active()

        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(self._doc2vec.docvecs.doctags),
                           redirect_stderr=True).start()
        count = 0
        for pk in self._doc2vec.docvecs.doctags:
            if not pk.startswith('j'):
                try:
                    paper = Paper.objects.get(pk=int(pk))
                    vector = self._doc2vec.docvecs[pk]
                    # normalize
                    vector /= np.linalg.norm(vector)
                    pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                               paper=paper)
                    pv.set_vector(vector)
                except Paper.DoesNotExist:
                    print('Paper {pk} did not match'.format(pk=pk))
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

    def save_all_vec_from_bulk(self):
        """Store inferred paper and journal vectors from training in db"""
        self.check_active()

        self.save_paper_vec_from_bulk()
        self.save_journal_vec_from_bulk()

    def infer_object(self, cls, pk, fields, **kwargs):

        # sanity checks
        self.check_active()
        # get object
        obj = cls.objects.get(id=pk)
        # tokenize fields
        obj_tokens = obj2tokens(obj, fields=fields)
        # Embed
        return self.embed(obj_tokens, **kwargs)

    def embed(self, tokens, alpha=0.1, min_alpha=0.001, passes=10):

        vector = self._doc2vec.infer_vector(tokens,
                                            alpha=alpha,
                                            min_alpha=min_alpha,
                                            steps=passes)
        # normalize
        vector /= np.linalg.norm(vector)

        return vector

    def activate(self):
        """Set model to active"""
        self.is_active = True
        self.save_db_only()

    def deactivate(self):
        """Set model to inactive"""
        self.is_active = False
        self.save_db_only()

    def check_active(self):
        if not self.is_active:
            raise ValidationError('model is not active')

    def get_words_vec(self, vec, topn=20):
        """retrieve closest top_n word from vector"""
        res = self._doc2vec.most_similar((vec,), topn=topn)
        closest_words = [r[0] for r in res]
        dist = [r[1] for r in res]
        return closest_words, dist

    def get_words_paper(self, paper_pk, topn=20):
        """retrieve closest top_n word from document vector"""
        pv = PaperVectors.objects.get(model=self, paper_id=paper_pk)
        vec = np.array(pv.get_vector())
        return self.get_words_vec(vec, topn=topn)

    def get_words_distribution(self, paper_pks, topn=10):
        """retrieve closest top_n words from all papers and build distribution"""
        words = []
        papers = Paper.objects.filter(pk__in=paper_pks)
        for paper in papers:
            vec = np.array(paper.vectors.get(model=self).get_vector())
            res = self._doc2vec.most_similar((vec,), topn=topn)
            words += [r[0] for r in res]

        # count occurences
        dist = collections.Counter(words)

        # remove single occurrence
        dist = dict((k, v) for k, v in dist.items() if v > 1)

        return dist

    def tasks(self, task_name, *args, **kwargs):
        """Dispatcher for Model tasks

        Dispatches tasks that can be run from Model. Task class is
        instantiated when celery worker starts up. Model instance is kept
        in-memory to avoid loading over-head.

        Args:
            task_name (string): Name of the method to call.
        """
        methods = [method for method in dir(self.__class__)
                   if callable(getattr(self.__class__,  method))]
        assert task_name in methods

        method = getattr(self, task_name)
        return method(*args, **kwargs)


class TextField(TimeStampedModel):
    """Store which fields of paper and related data are used in model training
    """

    text_field = models.CharField(max_length=100, choices=FIELDS_FOR_MODEL,
                                  null=False, blank=False)

    def __str__(self):
        return self.text_field


class PaperEngineManager(models.Manager):
    def load(self, **kwargs):
        """Get PaperEngine instance and load corresponding data"""
        obj = super(PaperEngineManager, self).get(**kwargs)
        try:
            pe_pkl_file = os.path.join(obj.PATH,
                                       '{0}.pkl'.format(obj.name))
            # if not on volume try download from s3
            if not os.path.isfile(pe_pkl_file):
                if getattr(settings, 'NLP_PE_BUCKET_NAME', ''):
                    obj.pull_from_s3()

            obj.data = joblib.load(pe_pkl_file)
        except EnvironmentError:  # OSError or IOError...
            raise

        return obj


class PaperEngine(PaperEngineScoringMixin, S3Mixin, TimeStampedModel):
    """PaperEngine class

    Store useful data and perform tasks for:
     - cosine similarity search
     - scoring

    Object is stored on S3 bucket.

    Data are dynamically updated and roll-out over DWELL_TIME
    """

    DWELL_TIME = 365    # days

    # For S3 Mixin
    # --------------------
    BUCKET_NAME = getattr(settings, 'NLP_PE_BUCKET_NAME', '')
    PATH = getattr(settings, 'NLP_PE_PATH', '')
    AWS_ACCESS_KEY_ID = getattr(settings, 'AWS_ACCESS_KEY_ID', '')
    AWS_SECRET_ACCESS_KEY = getattr(settings, 'AWS_SECRET_ACCESS_KEY', '')

    # Data structure
    # --------------------
    # Dictionary with following key, value pairs:
    # {
    #   'ids': [paper-id_1, ..., paper-id_N],
    #   'journal-ids': [journal-id_1, ..., journal-id_N],
    #   'authors-ids': [[auth-id_1_1, ... auth-id_1_P1],
    #                   [auth-id_2_1, ... auth-id_2_P2],
    #                   ...,
    #                   [auth-id_N_1, ... auth-id_N_PN]],
    #   'date':      [date_1, ..., date_N],
    #   'embedding': np.array([[v11, v12, ... v1Q],
    #                          ...,
    #                          [vN1, vN2, ... vNQ]])
    # }

    data = {'ids': [],
            'journal-ids': [],
            'authors-ids': [],
            'date': [],
            'embedding': np.empty(0),
            'altmetric': [],
            }

    # DB stored properties
    # --------------------
    model = models.ForeignKey(Model, related_name='pe')

    upload_state = models.CharField(max_length=3,
                                    choices=(('IDL', 'Idle'),
                                             ('ING', 'Uploading')),
                                    default='IDL ')
    # make instance active (use when linking task)
    is_active = models.BooleanField(default=False)

    # journal boost
    journal_boost = models.FloatField(default=0.05, null=True, blank=True)

    score_author_boost = models.FloatField(default=0.05)

    score_journal_boost = models.FloatField(default=0.05)

    score_altmetric_boost = models.FloatField(default=0.1)

    objects = PaperEngineManager()

    @property
    def name(self):
        return 'pe-{model_name}-id{id}'.format(
            model_name=self.model.name,
            id=self.id)

    @property
    def embedding_size(self):
        return self.model.size

    def __str__(self):
        return '{id}/{name}'.format(id=self.id, name=self.name)

    def activate(self):
        """Set model to active. Only on model can be active at once"""
        PaperEngine.objects.all() \
            .update(is_active=False)
        self.is_active = True
        self.save_db_only()

    def deactivate(self):
        """Set model to inactive"""
        self.is_active = False
        self.save_db_only()

    def save(self, *args, **kwargs):

        logger.info('Saving PE ({pk}/{name})'.format(pk=self.id,
                                                     name=self.name))

        # save to db (to get id)
        self.save_db_only()

        # save files to local volume
        if not os.path.exists(self.PATH):
            os.makedirs(self.PATH)
        pe_pkl_file = os.path.join(self.PATH,
                                   '{0}.pkl'.format(self.name))
        # use joblib for numpy array pickling optimization
        joblib.dump(self.data, pe_pkl_file)

        # push files to s3
        if self.BUCKET_NAME:
            self.upload_state = 'ING'
            self.save_db_only()
            self.push_to_s3(ext='pe')
            self.upload_state = 'IDL'

    def save_db_only(self, *args, **kwargs):
        super(PaperEngine, self).save(*args, **kwargs)

    def delete(self, s3=False, **kwargs):
        # delete on local volume
        try:
            # remove local
            rm_files = glob.glob(
                os.path.join(self.PATH, '{0}.*'.format(self.name)))
            for file in rm_files:
                os.remove(file)
        except IOError:
            pass

        # delete on S3
        if s3:
            try:
                if self.BUCKET_NAME:
                    self.delete_on_s3()
            except IOError:
                pass
        super(PaperEngine, self).delete(**kwargs)

    def full_update(self):
        """Full update of data structure"""

        logger.info('Updating PaperEngine ({pk}/{name}) - full update...'.format(
            pk=self.id, name=self.name))

        some_time_ago = (timezone.now() -
                         timezone.timedelta(days=self.DWELL_TIME)).date()

        self.data = {'ids': [],
                     'journal-ids': [],
                     'authors-ids': [],
                     'date': [],
                     'embedding': [],
                     'altmetric': [],
                     }

        # query on paper / vector
        q1 = Paper.objects.raw(
                "SELECT lp.id, "
                "		lp.journal_id, "
                "       LEAST(lp.date_ep, lp.date_pp, lp.date_fs) AS date_,"
                "		pv.vector,"
                "       aa.score "
                "FROM library_paper lp "
                "LEFT JOIN nlp_papervectors pv ON lp.id = pv.paper_id "
                "LEFT JOIN altmetric_altmetricmodel aa ON lp.id = aa.paper_id "
                "WHERE LEAST(date_ep, date_pp, date_fs) >= %s"
                "    AND pv.model_id=%s"
                "    AND lp.is_trusted=TRUE "
                "    AND lp.abstract <> ''"
                "    AND pv.vector IS NOT NULL "
                "ORDER BY date_ ASC", (some_time_ago, self.model.id)
                )

        for d in q1:
            self.data['ids'].append(d.id)
            self.data['journal-ids'].append(d.journal_id)
            self.data['date'].append(d.date_)
            self.data['altmetric'].append(d.score)
            self.data['embedding'].append(d.vector[:self.embedding_size])
        self.data['embedding'] = np.array(self.data['embedding'])

        # query on authors
        if self.data['ids']:
            values = ', '.join(['({0})'.format(i) for i in self.data['ids']])
            q2 = AuthorPaper.objects.raw(
                    "SELECT ap.id, "
                    "	    ap.paper_id, "
                    "		ap.author_id "
                    "FROM library_authorpaper ap "
                    "WHERE ap.paper_id IN (VALUES {0}) ".format(values))

            d2 = [(d.paper_id, d.author_id) for d in q2]
            dic = {}
            for k, v in d2:
                try:
                    dic[k].append(v)
                except KeyError:
                    dic[k] = [v]

            for i in self.data['ids']:
                if i in dic.keys():
                    self.data['authors-ids'].append(dic[i])
                else:
                    self.data['authors-ids'].append([])

        # Store
        self.save()

    def update(self):
        """Partial update of data structure with new entries only"""

        if not len(self.data['ids']):
            raise ValueError('Object is empty. run a full_update instead')

        logger.info('Updating PaperEngine ({pk}/{name}) - partial update...'.format(
            pk=self.id, name=self.name))

        some_time_ago = (timezone.now() -
                         timezone.timedelta(days=self.DWELL_TIME)).date()

        # query on paper / vector
        values = ', '.join(['({0})'.format(i) for i in self.data['ids']])
        q1 = Paper.objects.raw(
                "SELECT lp.id, "
                "		lp.journal_id, "
                "       LEAST(date_ep, date_pp, date_fs) AS date_,"
                "		pv.vector,"
                "       aa.score "
                "FROM library_paper lp "
                "LEFT JOIN nlp_papervectors pv ON lp.id = pv.paper_id "
                "LEFT JOIN altmetric_altmetricmodel aa ON lp.id = aa.paper_id "
                "WHERE LEAST(date_ep, date_pp, date_fs) >= %s"
                "    AND pv.model_id=%s"
                "    AND lp.is_trusted=TRUE "
                "    AND lp.abstract <> ''"
                "    AND pv.vector IS NOT NULL "
                "    AND lp.id NOT IN (VALUES {0}) "
                "ORDER BY date_ ASC".format(values), (some_time_ago, self.model.id,)
                )

        new_ids = []
        for d in q1:
            new_ids.append(d.id)
            self.data['ids'].append(d.id)
            self.data['journal-ids'].append(d.journal_id)
            self.data['date'].append(d.date_)
            self.data['altmetric'].append(d.score)
            np.vstack((self.data['embedding'], np.array(d.vector[:self.embedding_size])))

        # query on authors
        if new_ids:
            values = ', '.join(['({0})'.format(i) for i in new_ids])
            q2 = AuthorPaper.objects.raw(
                    "SELECT ap.id, "
                    "	    ap.paper_id, "
                    "		ap.author_id "
                    "FROM library_authorpaper ap "
                    "WHERE ap.paper_id IN (VALUES {0}) ".format(values))

            d2 = [(d.paper_id, d.author_id) for d in q2]
            dic = {}
            for k, v in d2:
                try:
                    dic[k].append(v)
                except KeyError:
                    dic[k] = [v]

            for i in new_ids:
                if i in dic.keys():
                    self.data['authors-ids'].append(dic[i])
                else:
                    self.data['authors-ids'].append([])

        # save
        self.save()

    def populate_neighbors(self, paper_id, time_lapse=-1):
        """Populate neighbors of paper"""

        neighbors_ids = self.get_knn(paper_id,
                                     time_lapse=time_lapse,
                                     k=settings.NLP_MAX_KNN_NEIGHBORS)

        pn, _ = PaperNeighbors.objects.get_or_create(pe=self,
                                                     time_lapse=time_lapse,
                                                     paper_id=paper_id)
        pn.set_neighbors(neighbors_ids)

        return neighbors_ids

    def get_clip_start(self, time_lapse):
        """Return index in data structure what corresponding to paper with
        date greater than now() - time_lapse

        Args:
            time_lapse (int): Number of days. if negative, returns 0

        """

        if time_lapse < 0:
            return 0
        else:
            cutoff_date = (timezone.now() -
                           timezone.timedelta(days=time_lapse)).date()
            if max(self.data['date']) > cutoff_date:
                return np.argmax(np.array(self.data['date']) > cutoff_date)
            else:
                return len(self.data['date'])

    def get_knn(self, paper_id, time_lapse=-1, k=1):

        pv = PaperVectors.objects \
                .filter(paper_id=paper_id, model=self.model) \
                .values('vector', 'paper__journal_id')[0]

        if paper_id in self.data['ids']:
            # increment k because knn_search will return input in the neighbors
            # set
            k += 1

        res_search = dict(self.knn_search(
            np.array(pv['vector'][:self.embedding_size]),
            time_lapse=time_lapse,
            top_n=k,
            journal_id=pv['paper__journal_id']))

        # Remove input from res_search
        try:
            res_search.pop(paper_id)
        except KeyError:
            pass

        return list(res_search.keys())

    def knn_search(self, seed, time_lapse=-1, top_n=5, journal_id=None):
        """"""

        # check seed
        if isinstance(seed, list):
            seed = np.array(seed).squeeze()
        if seed.ndim == 1:
            assert len(seed) == self.embedding_size
        else:
            assert seed.shape[1] == self.embedding_size

        # get clip
        clip_start = self.get_clip_start(time_lapse)

        # compute distance
        dists = np.dot(self.data['embedding'][clip_start:, :], seed)

        # Add journal bonus to distance
        if journal_id and self.journal_boost:
            jb = np.where(
                np.array(self.data['journal-ids'])[clip_start:] == journal_id,
                self.journal_boost,
                0)
            dists += jb

        best = matutils.argsort(dists, topn=top_n, reverse=True)

        # return [(paper_id, distances), ...]
        result = [(self.data['ids'][ind], float(dists[ind])) for ind in best]

        return result

    def get_partition(self, paper_ids, time_lapse=-1, k=1):
        """Return the k ids of neighbor papers of paper_ids"""

        # get seed data
        pvs = PaperVectors.objects \
                .filter(paper_id__in=paper_ids, model=self.model) \
                .values('vector', 'paper__journal_id')

        # Reshape
        mat = []
        journal_ids = []
        for pv in pvs:
            mat.append(pv['vector'][:self.embedding_size])
            journal_ids.append(pv['paper__journal_id'])
        mat = np.array(mat)

        # find clip
        clip_start = self.get_clip_start(time_lapse)

        # correct k for self-match
        k_correction = len([1 for id_ in paper_ids if id_ in self.data['ids']])
        k += k_correction

        # get partition
        res_search = self.partition_search(mat,
                                           clip_start=clip_start,
                                           top_n=k,
                                           journal_ids=journal_ids)
        # Remove inputs from res_search
        for id_ in paper_ids:
            res_search.remove(id_)

        return res_search

    def partition_search(self, seeds, clip_start=0, top_n=5,
                         clip_start_reverse=False,
                         journal_ids=None):
        """Return the unsorted list of paper ids that are in the top_n neighbors
        of vector defined as columns of 2d array seeds

        Args:
            seeds (np.array): 2d array, vector size x # of matches
        """

        # check seeds shape / type
        if isinstance(seeds, list) and any(isinstance(i, list) for i in seeds):
            seeds = np.array(seeds)
        assert seeds.shape[1] == self.embedding_size

        # Reverse clip_start if flagged
        if clip_start_reverse:
            clip_start = len(self.data['ids']) - clip_start

        # compute distance
        dists = np.dot(self.data['embedding'][clip_start:, :], seeds.T)

        # Add journal bonus to distance
        jbs = np.zeros(dists.shape)
        if journal_ids and self.journal_boost:
            for i, jid in enumerate(journal_ids):
                jbs[:, i] = np.where(
                    np.array(self.data['journal-ids'])[clip_start:] == jid,
                    self.journal_boost,
                    0)
        dists += jbs

        # reverse order
        dists = - dists
        # get partition
        arg_part = np.argpartition(dists, top_n + 1, axis=0)[:top_n:, :]

        result = [self.data['ids'][ind] for ind in arg_part.flatten()]
        return list(set(result))

    def get_recent_pks(self, time_lapse=30):
        """Return last found/published paper pk"""
        clip_start = self.get_clip_start(time_lapse)
        return self.data['ids'][clip_start:], self.data['ids'][clip_start:]

    def tasks(self, task_name, *args, **kwargs):
        """Dispatcher for PaperEngine tasks

        Dispatches tasks that can be run from PaperEngine. Task class is
        instantiated when celery worker starts up. PaperEnging data is kept
        in-memory to avoid loading over-head.

        Args:
            task_name (string): Name of the method to call.
        """
        methods = [method for method in dir(self.__class__)
                   if callable(getattr(self.__class__,  method))]
        assert task_name in methods

        method = getattr(self, task_name)
        return method(*args, **kwargs)


class ThreadEngineManager(models.Manager):

    def load(self, **kwargs):
        """Get ThreadEngine instance and load corresponding data"""
        obj = super(ThreadEngineManager, self).get(**kwargs)
        try:
            te_pkl_file = os.path.join(obj.PATH,
                                       '{0}.pkl'.format(obj.name))
            # if not on volume try download from s3
            if not os.path.isfile(te_pkl_file):
                if getattr(settings, 'NLP_TE_BUCKET_NAME', ''):
                    obj.pull_from_s3()

            obj.data = joblib.load(te_pkl_file)
        except EnvironmentError:  # OSError or IOError...
            raise

        return obj


class ThreadEngine(TimeStampedModel, S3Mixin):
    """ThreadEngine class

    Store useful data and perform tasks for:
     - cosine similarity search
     - scoring

    Object is stored on S3 bucket.

    Data are dynamically updated and roll-out over DWELL_TIME
    """

    DWELL_TIME = 365 * 5    # days

    # For S3 Mixin
    # --------------------
    BUCKET_NAME = getattr(settings, 'NLP_TE_BUCKET_NAME', '')
    PATH = getattr(settings, 'NLP_TE_PATH', '')
    AWS_ACCESS_KEY_ID = getattr(settings, 'AWS_ACCESS_KEY_ID', '')
    AWS_SECRET_ACCESS_KEY = getattr(settings, 'AWS_SECRET_ACCESS_KEY', '')

    # Data structure
    # --------------------
    # Dictionary with following key, value pairs:
    # {
    #   'ids': [thread-id_1, ..., thread-id_N],
    #   'paper-ids': [thread-id_1, ..., thread-id_N],
    #   'users-ids': [[user-id_1_1, ... user-id_1_P1],
    #                 [user-id_2_1, ... user-id_2_P2],
    #                   ...,
    #                 [user-id_N_1, ... user-id_N_PN]
    #                ],
    #   'date':      [date_1, ..., date_N],
    #   'embedding': np.array([[v11, v12, ... v1Q],
    #                          ...,
    #                          [vN1, vN2, ... vNQ]])
    # }
    #
    # Comment:
    #      - users-ids: first user is thread owner, others are users that posted


    data = {'ids': [],
            'paper-ids': [],
            'users-ids': [],
            'date': [],
            'embedding': np.empty(0)
            }

    # DB stored properties
    # --------------------
    model = models.ForeignKey(Model, related_name='te')

    upload_state = models.CharField(max_length=3,
                                    choices=(('IDL', 'Idle'),
                                             ('ING', 'Uploading')),
                                    default='IDL ')
    # make instance active (use when linking task)
    is_active = models.BooleanField(default=False)

    # user boost
    user_boost = models.FloatField(default=0.05, null=True, blank=True)

    objects = ThreadEngineManager()

    @property
    def name(self):
        return 'te-{model_name}-id{id}'.format(
            model_name=self.model.name,
            id=self.id)

    @property
    def embedding_size(self):
        return self.model.size

    def __str__(self):
        return '{id}/{name}'.format(id=self.id, name=self.name)

    def activate(self):
        """Set model to active. Only one model can be active at once"""
        ThreadEngine.objects.all() \
            .update(is_active=False)
        self.is_active = True
        self.save_db_only()

    def deactivate(self):
        """Set model to inactive"""
        self.is_active = False
        self.save_db_only()

    def save(self, *args, **kwargs):

        logger.info('Saving TE ({pk}/{name})'.format(pk=self.id,
                                                     name=self.name))

        # save to db (to get id)
        self.save_db_only()

        # save files to local volume
        if not os.path.exists(self.PATH):
            os.makedirs(self.PATH)
        te_pkl_file = os.path.join(self.PATH,
                                   '{0}.pkl'.format(self.name))
        # use joblib for numpy array pickling optimization
        joblib.dump(self.data, te_pkl_file)

        # push files to s3
        if self.BUCKET_NAME:
            self.upload_state = 'ING'
            self.save_db_only()
            self.push_to_s3(ext='te')
            self.upload_state = 'IDL'

    def save_db_only(self, *args, **kwargs):
        super(ThreadEngine, self).save(*args, **kwargs)

    def delete(self, s3=False, **kwargs):
        # delete on local volume
        try:
            # remove local
            rm_files = glob.glob(
                os.path.join(self.PATH, '{0}.*'.format(self.name)))
            for file in rm_files:
                os.remove(file)
        except IOError:
            pass

        # delete on S3
        if s3:
            try:
                if self.BUCKET_NAME:
                    self.delete_on_s3()
            except IOError:
                pass
        super(ThreadEngine, self).delete(**kwargs)

    def update(self):
        """Update of data structure"""

        logger.info('Updating ThreadEngine ({pk}/{name}) - full update...'.format(
            pk=self.id, name=self.name))

        some_time_ago = (timezone.now() -
                         timezone.timedelta(days=self.DWELL_TIME)).date()

        self.data = {'ids': [],
                     'paper-ids': [],
                     'users-ids': [],
                     'date': [],
                     'embedding': []}

        # query on thread / vector
        q1 = Thread.objects.raw(
                "SELECT tt.id,"
                "       tt.user_id,"
                "       tt.published_at,"
                "       tt.paper_id,"
                "		tv.vector "
                "FROM threads_thread tt "
                "LEFT JOIN nlp_threadvectors as tv ON tt.id = tv.thread_id "
                "WHERE tt.published_at >= %s"
                "    AND tv.model_id = %s"
                "    AND tv.vector IS NOT NULL "
                "ORDER BY tt.published_at ASC", (some_time_ago, self.model.id)
                )

        for d in q1:
            self.data['ids'].append(d.id)
            self.data['date'].append(d.published_at)
            self.data['paper-ids'].append(d.paper_id)
            self.data['users-ids'].append([d.user_id]) # add owner only for now
            self.data['embedding'].append(d.vector[:self.embedding_size])
        self.data['embedding'] = np.array(self.data['embedding'])

        # query on authors
        if self.data['ids']:
            values = ', '.join(['({0})'.format(i) for i in self.data['ids']])
            q2 = ThreadPost.objects.raw(
                    "SELECT tp.id, "
                    "	    tp.thread_id, "
                    "		tp.user_id "
                    "FROM threads_threadpost tp "
                    "WHERE tp.thread_id IN (VALUES {0}) ".format(values))

            d2 = [(d.thread_id, d.user_id) for d in q2]
            dic = {}
            for k, v in d2:
                try:
                    dic[k].append(v)
                except KeyError:
                    dic[k] = [v]

            for i, id_ in enumerate(self.data['ids']):
                if id_ in dic.keys():
                    self.data['users-ids'][i] += dic[id_]

        # Store
        self.save()

    def populate_neighbors(self, thread_id, time_lapse=-1):
        """Populate neighbors of thread"""

        neighbors_ids = self.get_knn(thread_id,
                                     time_lapse=time_lapse,
                                     k=settings.NLP_MAX_KNN_NEIGHBORS)

        pn, _ = ThreadNeighbors.objects.get_or_create(te=self,
                                                      time_lapse=time_lapse,
                                                      thread_id=thread_id)
        pn.set_neighbors(neighbors_ids)

        return neighbors_ids

    def get_clip_start(self, time_lapse):
        """Return index in data structure what corresponding to paper with
        date greater than now() - time_lapse

        Args:
            time_lapse (int): Number of days. if negative, returns 0

        """

        if time_lapse < 0:
            return 0
        else:
            cutoff_date = (timezone.now() -
                           timezone.timedelta(days=time_lapse)).date()
            if max(self.data['date']) > cutoff_date:
                return np.argmax(np.array(self.data['date']) > cutoff_date)
            else:
                return len(self.data['date'])

    def get_knn(self, thread_id, time_lapse=-1, k=1):

        tv = ThreadVectors.objects \
                .filter(thread_id=thread_id, model=self.model) \
                .values('vector')[0]

        if thread_id in self.data['ids']:
            # increment k because knn_search will return input in the neighbors
            # set
            k += 1

        res_search = dict(self.knn_search(
            np.array(tv['vector'][:self.embedding_size]),
            time_lapse=time_lapse,
            top_n=k))

        # Remove input from res_search
        try:
            res_search.pop(thread_id)
        except KeyError:
            pass

        return list(res_search.keys())

    def knn_search(self, seed, time_lapse=-1, top_n=5):
        """"""

        # check seed
        if isinstance(seed, list):
            seed = np.array(seed).squeeze()
        if seed.ndim == 1:
            assert len(seed) == self.embedding_size
        else:
            assert seed.shape[1] == self.embedding_size

        # get clip
        clip_start = self.get_clip_start(time_lapse)

        # compute distance
        dists = np.dot(self.data['embedding'][clip_start:, :], seed)

        best = matutils.argsort(dists, topn=top_n, reverse=True)

        # return [(thread_id, distances), ...]
        result = [(self.data['ids'][ind], float(dists[ind])) for ind in best]

        return result

    def get_partition(self, thread_ids, time_lapse=-1, k=1):
        """Return the k ids of neighbor threads of thread_id"""

        # get seed data
        pvs = ThreadVectors.objects \
                .filter(thread_id__in=thread_ids, model=self.model) \
                .values('vector')

        # Reshape
        mat = []
        for pv in pvs:
            mat.append(pv['vector'][:self.embedding_size])
        mat = np.array(mat)

        # find clip
        clip_start = self.get_clip_start(time_lapse)

        # correct k for self-match
        k_correction = len([1 for id_ in thread_ids if id_ in self.data['ids']])
        k += k_correction

        # get partition
        res_search = self.partition_search(mat,
                                           clip_start=clip_start,
                                           top_n=k)
        # Remove inputs from res_search
        for id_ in thread_ids:
            res_search.remove(id_)

        return res_search

    def partition_search(self, seeds, clip_start=0, top_n=5,
                         clip_start_reverse=False):
        """Return the unsorted list of thread ids that are in the top_n neighbors
        of vector defined as columns of 2d array seeds

        Args:
            seeds (np.array): 2d array, vector size x # of matches
        """

        # check seeds shape / type
        if isinstance(seeds, list) and any(isinstance(i, list) for i in seeds):
            seeds = np.array(seeds)
        assert seeds.shape[1] == self.embedding_size

        # Reverse clip_start if flagged
        if clip_start_reverse:
            clip_start = len(self.data['ids']) - clip_start

        # compute distance
        dists = np.dot(self.data['embedding'][clip_start:, :], seeds.T)

        # reverse order
        dists = - dists
        # get partition
        arg_part = np.argpartition(dists, top_n + 1, axis=0)[:top_n:, :]

        result = [self.data['ids'][ind] for ind in arg_part.flatten()]
        return list(set(result))

    def get_recent_pks(self, time_lapse=30):
        """Return last found/published paper pk"""
        clip_start = self.get_clip_start(time_lapse)
        return self.data['ids'][clip_start:], self.data['ids'][clip_start:]

    def tasks(self, task_name, *args, **kwargs):
        """Dispatcher for PaperEngine tasks

        Dispatches tasks that can be run from PaperEngine. Task class is
        instantiated when celery worker starts up. PaperEnging data is kept
        in-memory to avoid loading over-head.

        Args:
            task_name (string): Name of the method to call.
        """
        methods = [method for method in dir(self.__class__)
                   if callable(getattr(self.__class__,  method))]
        assert task_name in methods

        method = getattr(self, task_name)
        return method(*args, **kwargs)
