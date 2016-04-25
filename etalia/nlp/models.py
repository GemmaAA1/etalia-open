# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import logging
import glob
import pickle
from random import shuffle
import collections

import numpy as np
from gensim import matutils
from sklearn.externals import joblib
from progressbar import ProgressBar, Percentage, Bar, ETA
from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.conf import settings
from django.db.models import QuerySet
from django.utils import timezone
from django.core.validators import MaxValueValidator
from django.core.exceptions import ValidationError
from gensim.models import Doc2Vec

from gensim.models import Phrases

from etalia.core.models import TimeStampedModel
from etalia.core.constants import NLP_TIME_LAPSE_CHOICES
from etalia.core.utils import pad_vector, pad_neighbors
from etalia.library.models import Paper, Journal
from etalia.threads.models import Thread
from etalia.threads.constant import THREAD_TIME_LAPSE_CHOICES
from .constants import FIELDS_FOR_MODEL
from .utils import obj2tokens, TaggedDocumentsIterator, model_attr_getter, \
    model_attr_setter
from .mixins import S3Mixin
from .constants import NLP_JOURNAL_RATIO_CHOICES

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# THREADS/NLP RELATED
# ------------------------------------------------------------------------------

class ThreadVectors(TimeStampedModel):
    """ Table Thread - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None

    Use setter/getter functions to set and get vector
    """

    thread = models.ForeignKey(Thread,
                               related_name='vectors')

    model = models.ForeignKey('Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('thread', 'model')

    def __str__(self):
        return '{thread_pk}/{name}'.format(thread_pk=self.thread.pk,
                                           name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]


class ThreadNeighbors(TimeStampedModel):
    """Table for Thread neighbors"""
    # thread
    thread = models.ForeignKey(Thread)

    # Time lapse for neighbors match
    time_lapse = models.IntegerField(default=-1,
                                     choices=THREAD_TIME_LAPSE_CHOICES,
                                     verbose_name='Days from right now')

    # MostSimilarThread
    ms = models.ForeignKey('MostSimilarThread')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{ms}/{time_lapse}'.format(
            ms=self.ms.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.ms.model.size]

    class Meta:
        unique_together = ('time_lapse', 'thread', 'ms')


class ModelThreadMixin(object):
    """Mixin to Model to add Thread support"""
    THREAD_TOKENIZED_FIELDS = ['title', 'content']

    def infer_thread(self, thread_pk, **kwargs):
        """
        Infer vector for model and thread instance
        """
        vector = self.infer_object(Thread, thread_pk,
                                   self.THREAD_TOKENIZED_FIELDS, **kwargs)

        pv, _ = ThreadVectors.objects.get_or_create(model=self,
                                                    thread_id=thread_pk)

        # store
        pv.set_vector(vector.tolist())

        return thread_pk

    def infer_threads(self, thread_pks, **kwargs):
        """Infer vector for model and all thread in thread_pks list
        """
        # Check inputs
        if isinstance(thread_pks, models.QuerySet):
            thread_pks = list(thread_pks)
        if not isinstance(thread_pks, list):
            raise TypeError(
                'thread_pks must be list or QuerySet, found {0} instead'.format(
                    type(thread_pks)))

        # setup progressbar
        nb_pbar_updates = 100
        nb_threads = len(thread_pks)
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=nb_threads, redirect_stderr=True).start()

        for count, thread_pk in enumerate(thread_pks):
            self.infer_thread(thread_pk, **kwargs)
            if not (nb_threads // nb_pbar_updates) == 0:
                if not count % (nb_threads // nb_pbar_updates):
                    pbar.update(count)
        # close progress bar
        pbar.finish()

    def tasks(self, task, **kwargs):
        if task == 'infer_thread':
            thread_pk = kwargs.pop('thread_pk')
            return self.infer_thread(thread_pk, **kwargs)
        if task == 'infer_threads':
            thread_pks = kwargs.pop('thread_pks')
            return self.infer_threads(thread_pks, **kwargs)
        return super(ModelThreadMixin, self).tasks(task, **kwargs)

    def infer_object(self, cls, pk, fields, **kwargs):
        raise NotImplemented


# ------------------------------------------------------------------------------
# LIBRARY/NLP RELATED
# ------------------------------------------------------------------------------

class PaperVectors(TimeStampedModel):
    """ Table Paper - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None

    Use set_vector() to pad and set vector list
    """

    paper = models.ForeignKey(Paper, related_name='vectors')

    model = models.ForeignKey('Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('paper', 'model')

    def __str__(self):
        return '{paper_pk}/{name}'.format(paper_pk=self.paper.pk,
                                          name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]


class JournalVectors(TimeStampedModel):
    """ Table Journal - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None.

    Use set_vector() to pad and set vector list
    """

    journal = models.ForeignKey(Journal, related_name='vectors')

    model = models.ForeignKey('Model')

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    class Meta:
        unique_together = ('journal', 'model')

    def __str__(self):
        return '{short_title}/{name}' \
            .format(short_title=self.journal.short_title, name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]


class PaperNeighbors(TimeStampedModel):
    """ Table of matches nearest neighbors"""

    paper = models.ForeignKey(Paper, related_name='neighbors')

    ms = models.ForeignKey('MostSimilar')

    time_lapse = models.IntegerField(default=-1,
                                     choices=NLP_TIME_LAPSE_CHOICES,
                                     verbose_name='Days from right now')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{ms}/{time_lapse}'.format(
            ms=self.ms.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.ms.model.size]

    class Meta:
        unique_together = ('time_lapse', 'paper', 'ms')


class JournalNeighbors(TimeStampedModel):
    """ Table of matches nearest neighbors"""

    journal = models.ForeignKey(Journal, related_name='neighbors')

    ms = models.ForeignKey('MostSimilar')

    # Primary keys of the k-nearest neighbors matches
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{ms}/{time_lapse}'.format(
            ms=self.ms.name,
            time_lapse=self.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save()

    def get_neighbors(self):
        return self.neighbors[:self.ms.model.size]

    class Meta:
        unique_together = ('journal', 'ms')


class ModelLibraryMixin(object):
    """Mixin to Model to add Thread support"""

    PAPER_TOKENIZED_FIELDS = ['title', 'abstract']

    def infer_paper(self, paper_pk, **kwargs):
        """
        Infer vector for model and thread instance
        """
        vector = self.infer_object(Paper, paper_pk,
                                   self.PAPER_TOKENIZED_FIELDS, **kwargs)

        pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                   paper_id=paper_pk)

        # store
        pv.set_vector(vector.tolist())

        return paper_pk

    def infer_papers(self, paper_pks, **kwargs):
        """Infer vector for model and all thread in thread_pks list
        """
        # Check inputs
        if isinstance(paper_pks, models.QuerySet):
            paper_pks = list(paper_pks)
        if not isinstance(paper_pks, list):
            raise TypeError(
                'thread_pks must be list or QuerySet, found {0} instead'.format(
                    type(paper_pks)))

        # setup progressbar
        nb_pbar_updates = 100
        nb_papers = len(paper_pks)
        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=nb_papers, redirect_stderr=True).start()

        for count, paper_pk in enumerate(paper_pks):
            self.infer_paper(paper_pk, **kwargs)
            if not (nb_papers // nb_pbar_updates) == 0:
                if not count % (nb_papers // nb_pbar_updates):
                    pbar.update(count)
        # close progress bar
        pbar.finish()

    def tasks(self, task, **kwargs):
        if task == 'infer_paper':
            paper_pk = kwargs.pop('paper_pk')
            return self.infer_paper(paper_pk, **kwargs)
        if task == 'infer_papers':
            paper_pks = kwargs.pop('paper_pks')
            return self.infer_papers(paper_pks, **kwargs)
        return super(ModelLibraryMixin, self).tasks(task, **kwargs)

    def infer_object(self, cls, pk, fields, **kwargs):
        raise NotImplemented


# ------------------------------------------------------------------------------
# NLP MAIN
# ------------------------------------------------------------------------------
class ModelBase(TimeStampedModel):
    def tasks(self, task, **kwargs):
        """Task dispatcher"""
        if task == 'get_words_vec':
            vec = kwargs.pop('vec')
            return self.get_words_vec(vec, **kwargs)
        if task == 'get_words_paper':
            paper_pk = kwargs.pop('paper_pk')
            return self.get_words_paper(paper_pk, **kwargs)
        if task == 'infer_object':
            class_name = kwargs.pop('class_name')
            pk = kwargs.pop('pk')
            fields = kwargs.pop('fields')
            return self.infer_object(class_name, pk, fields, **kwargs)
        raise ValueError('Unknown task action: {0}'.format(task))

    class Meta:
        abstract = True


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
            ModelBase):
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
        self.build_most_similar()

    def build_most_similar(self):
        try:
            ms = MostSimilar.objects.get(model=self)
            ms.delete()
        except MostSimilar.DoesNotExist:
            pass
        ms = MostSimilar.objects.create(model=self)
        ms.full_update()

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


class TextField(TimeStampedModel):
    """Store which fields of paper and related data are used in model training
    """

    text_field = models.CharField(max_length=100, choices=FIELDS_FOR_MODEL,
                                  null=False, blank=False)

    def __str__(self):
        return self.text_field


class MostSimilarManager(models.Manager):
    def load(self, **kwargs):
        """Get MostSimilar instance and load corresponding data"""
        obj = super(MostSimilarManager, self).get(**kwargs)
        try:
            # if not on volume try download from s3
            if not os.path.isfile(os.path.join(settings.NLP_MS_PATH,
                                               '{0}.ms_data'.format(obj.name))):
                if getattr(settings, 'NLP_MS_BUCKET_NAME', ''):
                    obj.pull_from_s3()

            obj.data = joblib.load(os.path.join(settings.NLP_MS_PATH,
                                                '{0}.ms_data'.format(obj.name)))
            with open(
                    os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_ind2pk'),
                    'rb') as f:
                obj.index2pk = pickle.load(f)
            with open(os.path.join(settings.NLP_MS_PATH,
                                   obj.name + '.ms_ind2journalpk'), 'rb') as f:
                obj.index2journalpk = pickle.load(f)
            with open(os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_date'),
                      'rb') as f:
                obj.date = pickle.load(f)
        except EnvironmentError:  # OSError or IOError...
            raise

        return obj


class MostSimilar(TimeStampedModel, S3Mixin):
    """Most Similar class to perform k-nearest neighbors search based on
    simple dot product. This is working because paper vector are normalized"""

    # For S3 Mixin
    BUCKET_NAME = getattr(settings, 'NLP_MS_BUCKET_NAME', '')
    PATH = getattr(settings, 'NLP_MS_PATH', '')
    AWS_ACCESS_KEY_ID = getattr(settings, 'AWS_ACCESS_KEY_ID', '')
    AWS_SECRET_ACCESS_KEY = getattr(settings, 'AWS_SECRET_ACCESS_KEY', '')

    model = models.ForeignKey(Model, related_name='ms')

    upload_state = models.CharField(max_length=3,
                                    choices=(('IDL', 'Idle'),
                                             ('ING', 'Uploading')),
                                    default='IDL ')

    # 2D array of # matches x vector size
    data = np.empty(0)
    # data index to paper pk
    index2pk = []
    # index date
    date = []
    # data index to journal pk
    index2journalpk = []
    # journal ratio to weight vectors with
    journal_ratio = models.FloatField(default=0.0,
                                      choices=NLP_JOURNAL_RATIO_CHOICES)
    # make instance active (use when linking task)
    is_active = models.BooleanField(default=False)

    objects = MostSimilarManager()

    class Meta:
        unique_together = (('model', 'journal_ratio'),)

    @property
    def name(self):
        return '{model_name}-jr{journal_ratio:.0f}perc'.format(
            model_name=self.model.name,
            journal_ratio=self.journal_ratio * 100)

    def __str__(self):
        return '{id}/{name}'.format(id=self.id, name=self.name)

    def activate(self):
        """Set model to active. Only on model can be active at once"""
        MostSimilar.objects.all() \
            .update(is_active=False)
        self.is_active = True
        self.save_db_only()

    def deactivate(self):
        """Set model to inactive"""
        self.is_active = False
        self.save_db_only()

    def save(self, *args, **kwargs):

        logger.info('Saving MS ({pk}/{name})'.format(pk=self.id,
                                                     name=self.name))

        # save files to local volume
        if not os.path.exists(settings.NLP_MS_PATH):
            os.makedirs(settings.NLP_MS_PATH)
        # use joblib for numpy array pickling optimization
        joblib.dump(self.data,
                    os.path.join(settings.NLP_MS_PATH, self.name + '.ms_data'))
        with open(os.path.join(settings.NLP_MS_PATH, self.name + '.ms_ind2pk'),
                  'wb') as f:
            pickle.dump(self.index2pk, f)
        with open(os.path.join(settings.NLP_MS_PATH,
                               self.name + '.ms_ind2journalpk'), 'wb') as f:
            pickle.dump(self.index2journalpk, f)
        with open(os.path.join(settings.NLP_MS_PATH, self.name + '.ms_date'),
                  'wb') as f:
            pickle.dump(self.date, f)

        # push files to s3
        if self.BUCKET_NAME:
            self.upload_state = 'ING'
            self.save_db_only()
            self.push_to_s3(ext='ms')
            self.upload_state = 'IDL'
        # save to db
        self.save_db_only()

    def save_db_only(self, *args, **kwargs):
        super(MostSimilar, self).save(*args, **kwargs)

    def delete(self, s3=False, **kwargs):
        # delete on local volume
        try:
            # remove local
            rm_files = glob.glob(
                os.path.join(settings.NLP_MS_PATH, '{0}.ms*'.format(self.name)))
            for file in rm_files:
                os.remove(file)
        except IOError:
            pass

        # delete on amazon s3
        if s3:
            try:
                if self.BUCKET_NAME:
                    self.delete_on_s3()
            except IOError:
                pass
        super(MostSimilar, self).delete(**kwargs)

    def full_update(self):
        """Full update data for knn search
        """

        logger.info('Updating MS ({pk}/{name}) - fetching full data...'.format(
            pk=self.id, name=self.name))

        # size of the model space used for truncation of vector
        vec_size = self.model.size

        a_year_ago = (timezone.now() - timezone.timedelta(days=365)).date()

        # query
        data = list(PaperVectors.objects.raw(
            "SELECT nlp_papervectors.id, "
            "       nlp_papervectors.paper_id, "
            "       vector, "
            "       LEAST(date_ep, date_pp, date_fs) AS date,"
            "       journal_id "
            "        "
            "FROM nlp_papervectors LEFT JOIN library_paper "
            "ON nlp_papervectors.paper_id=library_paper.id "
            "WHERE LEAST(date_ep, date_pp, date_fs) >= %s"
            "    AND model_id=%s"
            "    AND is_trusted=TRUE "
            "    AND abstract <> ''"
            "ORDER BY date ASC", (a_year_ago, self.model.pk)))

        data_journal = dict(JournalVectors.objects \
                            .filter(model=self.model) \
                            .values_list('journal_id', 'vector'))

        # Reshape data
        nb_items = len(data)
        self.date = []
        self.index2pk = []
        self.index2journalpk = []
        self.data = np.zeros((nb_items, vec_size))
        for i, dat in enumerate(data):
            if dat.vector:
                self.date.append(dat.date)
                # store paper pk
                self.index2pk.append(dat.paper_id)
                # store journal pk
                self.index2journalpk.append(dat.journal_id)
                # build input matrix for fit
                if dat.journal_id in data_journal:
                    self.data[i, :] = (1 - self.journal_ratio) * np.array(
                        dat.vector[:vec_size]) + \
                                      self.journal_ratio * np.array(
                                          data_journal[dat.journal_id][
                                          :vec_size])
                else:
                    self.data[i, :] = np.array(dat.vector[:vec_size])

        # Store
        self.save()

    def update(self):
        """Update data for knn search for new paper
        """

        if not self.index2pk:
            raise ValueError('Looks like you should run full_update instead')

        logger.info('Updating MS ({pk}/{name}) - fetching new data...'.format(
            pk=self.id, name=self.name))

        # size of the model space used for truncation of vector
        vec_size = self.model.size

        a_year_ago = (timezone.now() - timezone.timedelta(days=365)).date()

        data = list(PaperVectors.objects.raw(
            "SELECT nlp_papervectors.id, "
            "       nlp_papervectors.paper_id, "
            "       vector, "
            "       LEAST(date_ep, date_pp, date_fs) AS date,"
            "       journal_id "
            "        "
            "FROM nlp_papervectors LEFT JOIN library_paper "
            "ON nlp_papervectors.paper_id=library_paper.id "
            "WHERE LEAST(date_ep, date_pp, date_fs) >= %s"
            "    AND model_id=%s"
            "    AND is_trusted=TRUE "
            "    AND abstract <> ''"
            "    AND nlp_papervectors.paper_id NOT IN %s"
            "ORDER BY date ASC",
            (a_year_ago, self.model.pk, tuple(self.index2pk))))

        data_journal = dict(JournalVectors.objects \
                            .filter(model=self.model) \
                            .values_list('journal_id', 'vector'))

        # Reshape data
        nb_items = len(data)
        date = []
        index2pk = []
        index2journalpk = []
        data_mat = np.zeros((nb_items, vec_size))
        for i, dat in enumerate(data):
            if dat.vector:
                # get min date
                date.append(dat.date)
                # store paper pk
                index2pk.append(dat.paper_id)
                # store journal pk
                index2journalpk.append(dat.journal_id)
                # build input matrix for fit
                if dat.journal_id in data_journal:
                    data_mat[i, :] = (1 - self.journal_ratio) * np.array(
                        dat.vector[:vec_size]) + \
                                     self.journal_ratio * np.array(
                                         data_journal[dat.journal_id][
                                         :vec_size])
                else:
                    data_mat[i, :] = np.array(dat.vector[:vec_size])

        # concatenate
        self.date += date
        self.index2pk += index2pk
        self.index2journalpk += index2journalpk
        self.data = np.vstack((self.data, data_mat))

        # reorder data by date
        idx = sorted(range(len(self.date)), key=self.date.__getitem__)
        self.date = [self.date[i] for i in idx]
        self.index2pk = [self.index2pk[i] for i in idx]
        self.index2journalpk = [self.index2journalpk[i] for i in idx]
        self.data = self.data[idx, :]

        # save
        self.save()

    def populate_neighbors(self, paper_pk, time_lapse=-1):
        """Populate neighbors of paper
        """

        neighbors_pks = self.get_knn(paper_pk, time_lapse=time_lapse,
                                     k=settings.NLP_MAX_KNN_NEIGHBORS)

        pn, _ = PaperNeighbors.objects.get_or_create(ms=self,
                                                     time_lapse=time_lapse,
                                                     paper_id=paper_pk)
        pn.set_neighbors(neighbors_pks)

        return neighbors_pks

    def get_clip_start(self, time_lapse):
        if time_lapse == -1:
            return 0
        else:
            cutoff_date = (
                timezone.now() - timezone.timedelta(days=time_lapse)).date()
            if max(np.array(self.date) > cutoff_date):  # at least one in True
                return np.argmax(np.array(self.date) > cutoff_date)
            else:
                return len(self.date)

    def get_knn(self, paper_pk, time_lapse=-1, k=1):

        pv = PaperVectors.objects \
            .filter(paper_id=paper_pk, model=self.model) \
            .select_related('paper__journal')[0]
        # get related journal vector
        jv = None
        if pv.paper.journal:
            try:
                jv = JournalVectors.objects \
                    .get(model=self.model, journal_id=pv.paper.journal.id)
            except JournalVectors.DoesNotExist:
                pass

        # build weighted vector
        if jv:
            vec = (1. - self.journal_ratio) * np.array(
                pv.vector[:self.model.size]) + \
                  self.journal_ratio * np.array(jv.vector[:self.model.size])
        else:
            vec = np.array(pv.vector[:self.model.size])

        res_search = self.knn_search(vec, time_lapse=time_lapse, top_n=k)

        neighbors_pks = [item[0] for item in res_search
                         if item[0] not in [paper_pk]]

        return neighbors_pks[:k]

    def knn_search(self, seed, time_lapse=-1, top_n=5):

        # check seed
        if isinstance(seed, list):
            seed = np.array(seed).squeeze()
        if seed.ndim == 1:
            assert len(seed) == self.model.size
        else:
            assert seed.shape[1] == self.model.size

        # get clip
        clip_start = self.get_clip_start(time_lapse)

        # compute distance
        dists = np.dot(self.data[clip_start:, :], seed)
        # clip index
        index2pk = self.index2pk[clip_start:]
        # sort (NB: +1 because likely will return input seed as closest item)
        best = matutils.argsort(dists, topn=top_n + 1, reverse=True)
        # return paper pk and distances
        result = [(index2pk[ind], float(dists[ind])) for ind in best]
        return result

    def get_partition(self, paper_pks, time_lapse=-1, k=1, journal_pks=None):
        """
        """

        # get seed data
        data = list(PaperVectors.objects \
                    .filter(paper_id__in=paper_pks, model=self.model) \
                    .select_related('paper__journal') \
                    .values('vector', 'paper__journal_id'))

        # Get corresponding journal data
        j_pks = [d.get('paper__journal_id') for d in data]
        data_journal = dict(JournalVectors.objects \
                            .filter(model=self.model, journal_id__in=j_pks) \
                            .values_list('journal_id', 'vector'))

        # Build the corresponding matrix data (paper weighted by journal in
        # agreement with the journal ratio current used in the active
        # MostSimilar object
        mat = np.zeros((self.model.size, len(data)))
        for i, d in enumerate(data):
            if 'paper__journal_id' in d and d[
                'paper__journal_id'] in data_journal:
                mat[:, i] = (1. - self.journal_ratio) * np.array(
                    d['vector'][:self.model.size]) + \
                            self.journal_ratio * np.array(
                                data_journal[d['paper__journal_id']][
                                :self.model.size])
            else:
                mat[:, i] = np.array(d['vector'][:self.model.size])

        # find clip
        clip_start = self.get_clip_start(time_lapse)

        # get partition
        res_search = self.partition_search(mat, clip_start=clip_start, top_n=k,
                                           journal_pks=journal_pks)

        # ignore paper_pks
        neighbors_pks_multi = [nei for nei in res_search if
                               nei not in paper_pks]
        return neighbors_pks_multi

    def partition_search(self, seeds, clip_start=0, top_n=5,
                         clip_start_reverse=False,
                         journal_pks=None):
        """Return the unsorted list of paper pk that are in the top_n neighbors
        of vector defined as columns of 2d array seeds

        Args:
            seeds (np.array): 2d array, vector size x # of matches
        """

        # check seeds
        if isinstance(seeds, list) and any(isinstance(i, list) for i in seeds):
            seeds = np.array(seeds)
        assert seeds.shape[0] == self.model.size

        # Reverse clip_start if flagged
        if clip_start_reverse:
            clip_start = len(self.index2pk) - clip_start

        # clip
        data = self.data[clip_start:, :]
        index2pk = self.index2pk[clip_start:]
        index2journalpk = self.index2journalpk[clip_start:]

        # filter by journal if defined
        if journal_pks:
            idx = list(map(lambda x: x in journal_pks, index2journalpk))
            data = data[np.array(idx), :]
            index2pk = [x for i, x in enumerate(index2pk) if idx[i]]
            # index2journalpk = [x for i, x in enumerate(index2journalpk) if idx[i]]

        # compute distance
        dists = np.dot(data, seeds)

        # reverse order
        dists = - dists
        # get partition
        arg_part = np.argpartition(dists, top_n + 1, axis=0)[:top_n:, :]
        # reshape and return unique
        result = [index2pk[ind] for ind in arg_part.flatten()]
        return list(set(result))

    def get_recent_pks(self, time_lapse=30):
        """Return last found/published paper pk"""
        clip_start = self.get_clip_start(time_lapse)
        return self.index2pk[clip_start:], self.date[clip_start:]

    def tasks(self, task, **kwargs):
        """Wrapper around MostSimilar tasks

        Use from tasks.py for calling task while object remains in-memory.
        Possibly an awkward design.

        This method dispatches tasks that can be run from MostSimilar.
        It is so because we want to have the MostSimilar instance in-memory
        to avoid over-head of loading from file. Given the Task class of Celery
        the work-around to define this wrapper that dispatches to sub tasks.

        Args:
            task (string): A string defining the task. 'update', 'populate_neighbors',
                'get_knn', 'knn_search', 'init'

            Optional kwargs:
            'populate_neighbors':
                paper_pk (int): The primary key of a Paper instance
                time_lapse (int): Days in the past, in NLP_TIME_LAPSE_CHOICES
            'get_knn':
                paper_pk (int): The primary key of a Paper instance
                time_lapse (int): Days in the past, in NLP_TIME_LAPSE_CHOICES
                k (int): number of neighbors
            'get_partition':
                paper_pks (list): List of primary key of Paper instances
                time_lapse (int): Days in the past, in NLP_TIME_LAPSE_CHOICES
                k (int): number of neighbors
            'knn_search':
                seed (np.array or list): An array of length corresponding to model
                clip_start (int): Where to start in self.data (see get_clip_start method)
                top_n (int): the top n results

        """

        if task == 'update':
            self.update()
        elif task == 'populate_neighbors':
            paper_pk = kwargs.pop('paper_pk')
            return self.populate_neighbors(paper_pk, **kwargs)
        elif task == 'get_knn':
            paper_pk = kwargs.pop('paper_pk')
            return self.get_knn(paper_pk, **kwargs)
        elif task == 'get_partition':
            paper_pks = kwargs.pop('paper_pks')
            return self.get_partition(paper_pks, **kwargs)
        elif task == 'knn_search':
            seed = kwargs.pop('seed')
            return self.knn_search(seed, **kwargs)
        elif task == 'get_recent_pks':
            return self.get_recent_pks(**kwargs)
        else:
            raise ValueError('Unknown task action: {0}'.format(task))


class MostSimilarThreadManager(models.Manager):
    def load(self, **kwargs):
        """Get MostSimilar instance and load corresponding data"""
        obj = super(MostSimilarThreadManager, self).get(**kwargs)
        try:
            # if not on volume try download from s3
            if not os.path.isfile(os.path.join(settings.NLP_MS_PATH,
                                               '{0}.ms_data'.format(obj.name))):
                if getattr(settings, 'NLP_MS_BUCKET_NAME', ''):
                    obj.pull_from_s3()

            obj.data = joblib.load(os.path.join(settings.NLP_MS_PATH,
                                                '{0}.ms_data'.format(obj.name)))
            with open(
                    os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_ind2pk'),
                    'rb') as f:
                obj.index2pk = pickle.load(f)
            with open(os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_date'),
                      'rb') as f:
                obj.date = pickle.load(f)
        except EnvironmentError:  # OSError or IOError...
            raise

        return obj


class MostSimilarThread(TimeStampedModel, S3Mixin):
    """Most Similar class to perform k-nearest neighbors search based on
    simple dot product. This is working because thread vector are normalized"""

    # For S3 Mixin
    BUCKET_NAME = getattr(settings, 'NLP_MS_BUCKET_NAME', '')
    PATH = getattr(settings, 'NLP_MS_PATH', '')
    AWS_ACCESS_KEY_ID = getattr(settings, 'AWS_ACCESS_KEY_ID', '')
    AWS_SECRET_ACCESS_KEY = getattr(settings, 'AWS_SECRET_ACCESS_KEY', '')

    model = models.ForeignKey(Model, related_name='ms_thread')

    upload_state = models.CharField(max_length=3,
                                    choices=(('IDL', 'Idle'),
                                             ('ING', 'Uploading')),
                                    default='IDL ')

    # 2D array of # matches x vector size
    data = np.empty(0)
    # data index to thread pk
    index2pk = []
    # index date
    date = []
    # make instance active (use when linking task)
    is_active = models.BooleanField(default=False)

    objects = MostSimilarThreadManager()

    @property
    def name(self):
        return 'thread-{model_name}'.format(
            model_name=self.model.name)

    def __str__(self):
        return '{id}/{name}'.format(id=self.id, name=self.name)

    def activate(self):
        """Set model to active. Only on model can be active at once"""
        MostSimilarThread.objects.all() \
            .update(is_active=False)
        self.is_active = True
        self.save_db_only()

    def deactivate(self):
        """Set model to inactive"""
        self.is_active = False
        self.save_db_only()

    def save(self, *args, **kwargs):

        logger.info('Saving MS Thread ({pk}/{name})'.format(pk=self.id,
                                                            name=self.name))

        # save files to local volume
        if not os.path.exists(settings.NLP_MS_PATH):
            os.makedirs(settings.NLP_MS_PATH)
        # use joblib for numpy array pickling optimization
        joblib.dump(self.data,
                    os.path.join(settings.NLP_MS_PATH, self.name + '.ms_data'))
        with open(os.path.join(settings.NLP_MS_PATH, self.name + '.ms_ind2pk'),
                  'wb') as f:
            pickle.dump(self.index2pk, f)
        with open(os.path.join(settings.NLP_MS_PATH, self.name + '.ms_date'),
                  'wb') as f:
            pickle.dump(self.date, f)

        # push files to s3
        if self.BUCKET_NAME:
            self.upload_state = 'ING'
            self.save_db_only()
            self.push_to_s3(ext='ms')
            self.upload_state = 'IDL'
        # save to db
        self.save_db_only()

    def save_db_only(self, *args, **kwargs):
        super(MostSimilarThread, self).save(*args, **kwargs)

    def delete(self, s3=False, **kwargs):
        # delete on local volume
        try:
            # remove local
            rm_files = glob.glob(
                os.path.join(settings.NLP_MS_PATH, '{0}.ms*'.format(self.name)))
            for file in rm_files:
                os.remove(file)
        except IOError:
            pass

        # delete on amazon s3
        if s3:
            try:
                if self.BUCKET_NAME:
                    self.delete_on_s3()
            except IOError:
                pass
        super(MostSimilarThread, self).delete(**kwargs)

    def full_update(self):
        """Full update data for knn search
        """

        logger.info('Updating MS ({pk}/{name}) - fetching full data...'.format(
            pk=self.id, name=self.name))

        # size of the model space used for truncation of vector
        vec_size = self.model.size

        # query
        data = list(ThreadVectors.objects.raw(
            "SELECT nlp_threadvectors.id, "
            "       nlp_threadvectors.thread_id, "
            "       vector, "
            "       published_at"
            "        "
            "FROM nlp_threadvectors LEFT JOIN threads_thread "
            "ON nlp_threadvectors.thread_id=threads_thread.id "
            "WHERE model_id=%s "
            "    AND published_at IS NOT NULL "
            "ORDER BY published_at ASC", (self.model.pk, )))

        # Reshape data
        nb_items = len(data)
        self.date = []
        self.index2pk = []
        self.data = np.zeros((nb_items, vec_size))
        for i, dat in enumerate(data):
            if dat.vector:
                self.date.append(dat.published_at)
                # store thread pk
                self.index2pk.append(dat.thread_id)
                # build input matrix for fit
                self.data[i, :] = np.array(dat.vector[:vec_size])

        # Store
        self.save()

    def update(self):
        """Update data for knn search for new thread
        """

        if not self.index2pk:
            raise ValueError('Looks like you should run full_update instead')

        logger.info('Updating MS ({pk}/{name}) - fetching new data...'.format(
            pk=self.id, name=self.name))

        # size of the model space used for truncation of vector
        vec_size = self.model.size

        data = list(ThreadVectors.objects.raw(
            "SELECT nlp_threadvectors.id, "
            "       nlp_threadvectors.thread_id, "
            "       vector, "
            "       published_at"
            "        "
            "FROM nlp_threadvectors LEFT JOIN threads_thread "
            "ON nlp_threadvectors.thread_id=threads_thread.id "
            "WHERE model_id=%s"
            "    AND published_at IS NOT NULL "
            "    AND nlp_threadvectors.thread_id NOT IN %s"
            "ORDER BY published_at ASC",
            (self.model.pk, tuple(self.index2pk))))

        # Reshape data
        nb_items = len(data)
        date = []
        index2pk = []
        data_mat = np.zeros((nb_items, vec_size))
        for i, dat in enumerate(data):
            if dat.vector:
                # get min date
                date.append(dat.published_at)
                # store thread pk
                index2pk.append(dat.thread_id)
                # build input matrix for fit
                data_mat[i, :] = np.array(dat.vector[:vec_size])

        # concatenate
        self.date += date
        self.index2pk += index2pk
        self.data = np.vstack((self.data, data_mat))

        # reorder data by date
        idx = sorted(range(len(self.date)), key=self.date.__getitem__)
        self.date = [self.date[i] for i in idx]
        self.index2pk = [self.index2pk[i] for i in idx]
        self.data = self.data[idx, :]

        # save
        self.save()

    def populate_neighbors(self, thread_pk, time_lapse=-1):
        """Populate neighbors of thread
        """

        neighbors_pks = self.get_knn(thread_pk, time_lapse=time_lapse,
                                     k=settings.NLP_MAX_KNN_NEIGHBORS)

        pn, _ = ThreadNeighbors.objects.get_or_create(ms=self,
                                                      time_lapse=time_lapse,
                                                      thread_id=thread_pk)
        pn.set_neighbors(neighbors_pks)

        return neighbors_pks

    def get_clip_start(self, time_lapse):
        if time_lapse == -1:
            return 0
        else:
            cutoff_date = (
                timezone.now() - timezone.timedelta(days=time_lapse))
            if max(np.array(self.date) > cutoff_date):  # at least one in True
                return np.argmax(np.array(self.date) > cutoff_date)
            else:
                return len(self.date)

    def get_knn(self, thread_pk, time_lapse=-1, k=1):

        pv = ThreadVectors.objects \
            .filter(thread_id=thread_pk, model=self.model)[0]

        # build weighted vector
        vec = np.array(pv.vector[:self.model.size])

        res_search = self.knn_search(vec, time_lapse=time_lapse, top_n=k)

        neighbors_pks = [item[0] for item in res_search
                         if item[0] not in [thread_pk]]

        return neighbors_pks[:k]

    def knn_search(self, seed, time_lapse=-1, top_n=5):

        # check seed
        if isinstance(seed, list):
            seed = np.array(seed).squeeze()
        if seed.ndim == 1:
            assert len(seed) == self.model.size
        else:
            assert seed.shape[1] == self.model.size

        # get clip
        clip_start = self.get_clip_start(time_lapse)

        # compute distance
        dists = np.dot(self.data[clip_start:, :], seed)
        # clip index
        index2pk = self.index2pk[clip_start:]
        # sort (NB: +1 because likely will return input seed as closest item)
        best = matutils.argsort(dists, topn=top_n + 1, reverse=True)
        # return thread pk and distances
        result = [(index2pk[ind], float(dists[ind])) for ind in best]
        return result

    def get_partition(self, thread_pks, time_lapse=-1, k=1):
        """
        """

        # get seed data
        data = list(ThreadVectors.objects \
                    .filter(thread_id__in=thread_pks, model=self.model) \
                    .values('vector'))

        # Build the corresponding matrix data
        mat = np.zeros((self.model.size, len(data)))
        for i, d in enumerate(data):
                mat[:, i] = np.array(d['vector'][:self.model.size])

        # find clip
        clip_start = self.get_clip_start(time_lapse)

        # get partition
        res_search = self.partition_search(mat, clip_start=clip_start, top_n=k)

        # ignore thread_pks
        neighbors_pks_multi = [nei for nei in res_search if
                               nei not in thread_pks]
        return neighbors_pks_multi

    def partition_search(self, seeds, clip_start=0, top_n=5,
                         clip_start_reverse=False):
        """Return the unsorted list of thread pk that are in the top_n neighbors
        of vector defined as columns of 2d array seeds

        Args:
            seeds (np.array): 2d array, vector size x # of matches
        """

        # check seeds
        if isinstance(seeds, list) and any(isinstance(i, list) for i in seeds):
            seeds = np.array(seeds)
        assert seeds.shape[0] == self.model.size

        # Reverse clip_start if flagged
        if clip_start_reverse:
            clip_start = len(self.index2pk) - clip_start

        # clip
        data = self.data[clip_start:, :]
        index2pk = self.index2pk[clip_start:]

        # compute distance
        dists = np.dot(data, seeds)

        # reverse order
        dists = - dists
        # get partition
        arg_part = np.argpartition(dists, top_n + 1, axis=0)[:top_n:, :]
        # reshape and return unique
        result = [index2pk[ind] for ind in arg_part.flatten()]
        return list(set(result))

    def get_recent_pks(self, time_lapse=30):
        """Return last found/published thread pk"""
        clip_start = self.get_clip_start(time_lapse)
        return self.index2pk[clip_start:], self.date[clip_start:]

    def tasks(self, task, **kwargs):
        """Wrapper around MostSimilar tasks

        Use from tasks.py for calling task while object remains in-memory.
        Possibly an awkward design.

        This method dispatches tasks that can be run from MostSimilar.
        It is so because we want to have the MostSimilar instance in-memory
        to avoid over-head of loading from file. Given the Task class of Celery
        the work-around to define this wrapper that dispatches to sub tasks.

        Args:
            task (string): A string defining the task. 'update', 'populate_neighbors',
                'get_knn', 'knn_search', 'init'

            Optional kwargs:
            'populate_neighbors':
                thread_pk (int): The primary key of a Paper instance
                time_lapse (int): Days in the past, in NLP_TIME_LAPSE_CHOICES
            'get_knn':
                thread_pk (int): The primary key of a Paper instance
                time_lapse (int): Days in the past, in NLP_TIME_LAPSE_CHOICES
                k (int): number of neighbors
            'get_partition':
                thread_pks (list): List of primary key of Paper instances
                time_lapse (int): Days in the past, in NLP_TIME_LAPSE_CHOICES
                k (int): number of neighbors
            'knn_search':
                seed (np.array or list): An array of length corresponding to model
                clip_start (int): Where to start in self.data (see get_clip_start method)
                top_n (int): the top n results

        """

        if task == 'update':
            self.update()
        elif task == 'populate_neighbors':
            thread_pk = kwargs.pop('thread_pk')
            return self.populate_neighbors(thread_pk, **kwargs)
        elif task == 'get_knn':
            thread_pk = kwargs.pop('thread_pk')
            return self.get_knn(thread_pk, **kwargs)
        elif task == 'get_partition':
            thread_pks = kwargs.pop('thread_pks')
            return self.get_partition(thread_pks, **kwargs)
        elif task == 'knn_search':
            seed = kwargs.pop('seed')
            return self.knn_search(seed, **kwargs)
        elif task == 'get_recent_pks':
            return self.get_recent_pks(**kwargs)
        else:
            raise ValueError('Unknown task action: {0}'.format(task))
