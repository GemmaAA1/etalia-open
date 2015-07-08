import os
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.conf import settings
from core.models import TimeStampedModel

from gensim.models import Doc2Vec

from .constants import MODEL_STATES

from library.models import Paper
# Create your models here.

class ModelManager(models.Manager):
    def create(self, **kwargs):
        model = self.model(**kwargs)
        model.save(using=self._db)
        model.doc2vec = Doc2Vec(
            dm=model.dm,
            size=model.size,
            window=model.window,
            alpha=model.alpha,
            seed=model.seed,
            min_count=model.min_count,
            max_vocab_size=model.max_vocab_size,
            sample=model.sample,
            workers=model.workers,
            hs=model.hs,
            negative=model.negative,
            dm_mean=model.dm_mean,
            dm_concat=model.dm_concat,
            dm_tag_count=model.dm_tag_count,
            dbow_words=model.dbow_words
        )
        model.save()
        return model

    def get(self, **kwargs):
        model = super(ModelManager, self).get(**kwargs)
        try:
            model.doc2vec = Doc2Vec.load(os.path.join(model.doc2vec_path,
                                         '{0}.mod'.format(model.name)))
        except FileNotFoundError:
            raise FileNotFoundError
        return model


class Model(TimeStampedModel):

    name = models.CharField(max_length=128, blank=False, null=False,
                            unique=True)

    doc2vec_path = models.CharField(
        max_length=256,
        default=settings.NLP_DOC2VEC_PATH)

    status = models.CharField(max_length=3, choices=MODEL_STATES, default='UNT')

    objects = ModelManager()

    # Model parameters
    # ----------------------
    # `dm` defines the training algorithm. By default (`dm=1`), 'distributed
    # memory' (PV-DM) is used.
    # Otherwise, `distributed bag of words` (PV-DBOW) is employed.
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

    doc2vec = Doc2Vec()

    def __str__(self):
        return self.name

    def build_vocab(self, documents):
        self.update_status('VOC')
        self.doc2vec.build_vocab(documents)
        self.update_status('IDL')

    def train(self, documents):
        self.update_status('TRA')
        self.doc2vec.train(documents)
        self.update_status('IDL')

    def save_model_only(self, *args, **kwargs):
        super(Model, self).save(*args, **kwargs)

    def save(self, *args, **kwargs):
        self.update_status('SAV')
        super(Model, self).save(*args, **kwargs)
        self.doc2vec.save(os.path.join(self.doc2vec_path, '{0}.mod'.format(self.name)))
        self.update_status('IDL')

    def delete(self, *args, **kwargs):
        try:
            os.remove(os.path.join(self.doc2vec_path, '{0}.mod'.format(self.name)))
        except FileNotFoundError:
            pass
        super(Model, self).delete(*args, **kwargs)

    def update_status(self, state):
        self.status = state
        self.save_model_only()


class Fingerprint(TimeStampedModel):
    paper = models.ForeignKey(Paper)

    mod = models.ForeignKey(Model)

    vec = ArrayField(models.FloatField())
