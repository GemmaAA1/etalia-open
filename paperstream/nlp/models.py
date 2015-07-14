import os
from random import shuffle
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.conf import settings
from core.models import TimeStampedModel
from progressbar import ProgressBar, Percentage, Bar, ETA

from gensim.models import Doc2Vec, Phrases

from .constants import MODEL_STATES
from .exceptions import StatusError

from library.models import Paper, Journal
# Create your models here.

class ModelManager(models.Manager):
    def create(self, **kwargs):
        model = self.model(**kwargs)
        model.init_doc2vec()
        model.save(using=self._db)
        return model

    def get(self, **kwargs):
        model = super(ModelManager, self).get(**kwargs)
        try:
            model.doc2vec = Doc2Vec.load(os.path.join(model.doc2vec_path,
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
    >>> model = Model.objects.create(name='test') # with default parameters
    2) Prepare and dump data.
    >>> papers = Papers.objects.all()
    >>> dumper = DumpPaperData(to=model.data_path)
    >>> dumper.dump(papers)
    3) Build vocab (with a phraser to join common n-gram word. ie new_york)
    >>> docs_l = WordListIterator('nlp/data')
    >>> phraser = Phrases(docs_l)
    >>> docs = TaggedDocumentsIterator('nlp/data', phraser=phraser)
    >>> model.build_vocab(docs)
    4) Train, save and set_active
    >>> model.train(docs, iteration=1)
    >>> model.save()
    >>> model.set_active()
    5) Populate library is needed
    >>> model.populate_library()
    """

    name = models.CharField(max_length=128, blank=False, null=False,
                            unique=True)

    doc2vec_path = models.CharField(
        max_length=256,
        default=settings.NLP_DOC2VEC_PATH)

    is_active = models.BooleanField(default=False)

    # document2vector instance from gensim
    doc2vec = Doc2Vec()

    objects = ModelManager()

    # Model parameters
    # ----------------------
    # `dm` defines the training algorithm. By default (`dm=1`), 'distributed
    # memory' (PV-DM) is used. Otherwise, `distributed bag of words` (PV-DBOW) is employed.
    dm = models.IntegerField(default=1)

    # `size` is the dimensionality of the feature vectors.
    size = models.IntegerField(default=100)

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
    hs = models.IntegerField(default=0)

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

    def build_vocab(self, documents):
        self.doc2vec.build_vocab(documents)

    def train(self, documents, passes=10):
        # Using explicit multiple-pass, alpha-reduction approach as sketched in
        # gensim doc2vec blog post
        alpha, min_alpha, passes = (0.025, 0.001, passes)
        alpha_delta = (alpha - min_alpha) / passes
        for epoch in range(passes):
            print('Epoch {epoch}/{passes}\n'.format(epoch=epoch, passes=passes))
            shuffle(documents)
            self.alpha = alpha
            self.doc2vec.alpha = alpha
            self.doc2vec.min_alpha = alpha
            self.doc2vec.train(documents)
            alpha -= alpha_delta

    def save_db_only(self, *args, **kwargs):
        super(Model, self).save(*args, **kwargs)

    def save(self, *args, **kwargs):
        self.full_clean()
        super(Model, self).save(*args, **kwargs)
        if not os.path.exists(self.doc2vec_path):
            os.makedirs(self.doc2vec_path)
        self.doc2vec.save(
            os.path.join(self.doc2vec_path, '{0}.mod'.format(self.name)))

    def delete(self, *args, **kwargs):
        try:
            os.remove(
                os.path.join(self.doc2vec_path, '{0}.mod'.format(self.name)))
        except FileNotFoundError:
            pass
        super(Model, self).delete(*args, **kwargs)

    def clean(self):
        """Check concordance of model field attributes and doc2vec attributes
        """
        for attr in self.model_arguments:
            if not getattr(self, attr[0]) == getattr(self.doc2vec, attr[1]):
                raise ValueError('{attr} does not match'.format(attr=attr[0]))

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

    def populate_library(self):

        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(self.doc2vec.docvecs.doctags),
                           redirect_stderr=True).start()
        count = 0
        # TODO: replace with bulk update / create ?
        for pk in self.doc2vec.docvecs.doctags:

            if pk[0] == 'j':
                try:
                    journal = Journal.objects.get(pk=int(pk[2:]))
                    vector = self.doc2vec.docvecs[pk].tolist()
                    pv, _ = JournalVectors.objects.get_or_create(model=self,
                                                                 journal=journal)
                    pv.vector = vector
                    pv.save()
                except Journal.DoesNotExist:
                    print('Journal {pk} did not match'.format(pk=pk))
            else:
                try:
                    paper = Paper.objects.get(pk=int(pk))
                    vector = self.doc2vec.docvecs[pk].tolist()
                    pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                               paper=paper)
                    pv.vector = vector
                    pv.save()
                except Paper.DoesNotExist:
                    print('Paper {pk} did not match'.format(pk=pk))
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

    def set_active(self):
        self.is_active = True
        self.save_model_only()

    class Meta:
        ordering = ['name', ]


class PaperVectors(TimeStampedModel):

    paper = models.ForeignKey(Paper)

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True), null=True)

    class meta:
        unique_together = ('paper', 'model')

    def __str__(self):
        return '{pk} with {name}'.format(pk=self.paper.pk, name=self.model.name)


class JournalVectors(TimeStampedModel):

    journal = models.ForeignKey(Journal)

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True), null=True)

    class meta:
        unique_together = ('journal', 'model')

    def __str__(self):
        return '{pk} with {name}'.format(pk=self.journal.pk, name=self.model.name)

