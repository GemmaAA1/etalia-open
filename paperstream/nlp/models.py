# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import logging
import glob
from random import shuffle
import numpy as np
from timeit import time
import boto
from boto.s3.key import Key
import tarfile

from sklearn.externals import joblib

from django.db import models, transaction
from django.db.models import Q, QuerySet
from django.contrib.postgres.fields import ArrayField
from django.conf import settings
from django.utils import timezone
from django.core.validators import MaxValueValidator, MinValueValidator
from django.core.exceptions import ValidationError

from progressbar import ProgressBar, Percentage, Bar, ETA
from gensim.models import Doc2Vec
from gensim.models import Phrases

from .constants import MODEL_STATES, FIELDS_FOR_MODEL, LSHF_STATES
from .utils import paper2tokens, TaggedDocumentsIterator, MyLSHForest, \
    model_attr_getter, model_attr_setter
from .exceptions import InvalidState
from .mixins import S3Mixin

from paperstream.core.models import TimeStampedModel
from paperstream.core.utils import pad_vector, pad_neighbors
from paperstream.library.models import Paper, Journal
from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES


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
            raise ValidationError('<text_fields> not subset of FIELDS_FOR_MODEL')

        for text_field in text_fields:
            fuim, _ = TextField.objects.get_or_create(text_field=text_field)
            obj.text_fields.add(fuim)

        obj.full_clean()
        obj.save_db_only(using=self.db)

        # Init doc2vec from gensim
        obj.init_doc2vec()

        obj.save(using=self.db)

        return obj

    def get(self, **kwargs):
        print('Warning: Not loading Doc2Vec. Use load() instead')
        return super(ModelManager, self).get(**kwargs)

    def load(self, **kwargs):
        obj = super(ModelManager, self).get(**kwargs)
        try:
            # if not on volume try download from s3
            if not os.path.isfile(os.path.join(settings.NLP_MODELS_PATH,
                                               '{}.mod'.format(obj.name))):
                obj.download_from_s3()
            obj._doc2vec = Doc2Vec.load(os.path.join(
                settings.NLP_MODELS_PATH, '{0}.mod'.format(obj.name)))
        except EnvironmentError:      # OSError or IOError...
            raise
        return obj


class Model(TimeStampedModel, S3Mixin):
    """Natural Language Processing Class based on Doc2Vec from Gensim
    """

    # For S3 Mixin
    BUCKET_NAME = settings.NLP_MODELS_BUCKET_NAME
    PATH = settings.NLP_MODELS_PATH
    AWS_ACCESS_KEY_ID = settings.DJANGO_AWS_ACCESS_KEY_ID
    AWS_SECRET_ACCESS_KEY = settings.DJANGO_AWS_SECRET_ACCESS_KEY

    name = models.CharField(max_length=128, blank=False, null=False,
                            unique=True)

    # when trained, model becomes active
    is_active = models.BooleanField(default=False)

    # doc2vec instance from gensim
    _doc2vec = Doc2Vec()

    # fields
    text_fields = models.ManyToManyField('TextField',
                                         related_name='text_fields')

    objects = ModelManager()

    # Model parameters (mirroring Doc2Vec attributes)
    # ----------------------
    # `dm` defines the training algorithm. By default (`dm=1`), 'distributed
    # memory' (PV-DM) is used. Otherwise, `distributed bag of words` (PV-DBOW)
    # is employed.
    _dm = models.IntegerField(default=1)

    # `size` is the dimensionality of the feature vectors.
    _size = models.IntegerField(default=128,
            validators=[MaxValueValidator(settings.NLP_MAX_VECTOR_SIZE)])

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
        super(Model, self).save(*args, **kwargs)
        if not os.path.exists(settings.NLP_MODELS_PATH):
            os.makedirs(settings.NLP_MODELS_PATH)
        self._doc2vec.save(
            os.path.join(settings.NLP_MODELS_PATH,
                         '{0}.mod'.format(self.name)))

    def delete(self, *args, **kwargs):
        try:
            for filename in glob.glob(os.path.join(settings.NLP_MODELS_PATH,
                                                   '{0}.mod*'.format(self.name))):
                os.remove(filename)
        except IOError:
            pass
        super(Model, self).delete(*args, **kwargs)

    def dump(self, papers, data_path=None):
        """Dump papers data to pre-process text files.
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
                           maxval=int(tot/sub_update_step),
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
            line_val = paper2tokens(paper, fields=text_fields)

            # write to file
            # fid.write(u' '.join(line_val).strip().encode('utf-8'))
            fid.write(u' '.join(line_val).strip())

            # write new line
            fid.write('\n')

            # update progress bar
            if not count % sub_update_step:
                pbar.update(count/sub_update_step)
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
            self.save()
            alpha -= alpha_delta

    def build_vocab_and_train(self, **kwargs):
        """Build vocabulary (with phraser is set) and train model"""
        # Build vocabulary and phraser
        docs = self.load_documents()
        phraser = self.build_phraser(docs, **kwargs)
        docs = self.load_documents(phraser=phraser)
        self.build_vocab(docs)
        # Train, save and activate
        self.train(docs, **kwargs)
        self.activate()

    def build_lshs(self):
        """Build Local Sensitive Hashing data structure"""
        LSH.objects.filter(model=self).delete()
        # Build time-lapse related LSH
        for time_lapse, _ in NLP_TIME_LAPSE_CHOICES:
            LSH.objects.create(model=self,
                               time_lapse=time_lapse)

    def propagate(self):
        """Save Journal and Paper vectors and build LSHs
        (To be run after model training)"""
        logging.info('{} Populate journals...'.format(self.name))
        self.save_journal_vec_from_bulk()
        logging.info('{} Populate papers...'.format(self.name))
        self.save_paper_vec_from_bulk()
        logging.info('{} Build LSH...'.format(self.name))
        self.build_lshs()

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
                    vector = self._doc2vec.docvecs[pk].tolist()
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
        self.check_active()

        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(self._doc2vec.docvecs.doctags),
                           redirect_stderr=True).start()
        count = 0
        for pk in self._doc2vec.docvecs.doctags:
            if not pk.startswith('j'):
                try:
                    paper = Paper.objects.get(pk=int(pk))
                    vector = self._doc2vec.docvecs[pk].tolist()
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
        """Store inferred paper and journal vectors from training in db"""
        self.check_active()

        self.save_paper_vec_from_bulk()
        self.save_journal_vec_from_bulk()

    def infer_paper(self, paper_pk, alpha=0.1, min_alpha=0.001, passes=5):  # seed=seed)
        """Infer model vector for paper"""

        self.check_active()

        text_fields = [tx.text_field for tx in self.text_fields.all()]

        pv, _ = PaperVectors.objects.get_or_create(model=self,
                                                   paper_id=paper_pk)
        paper = Paper.objects.get(id=paper_pk)
        doc_words = paper2tokens(paper, fields=text_fields)

        # NB: infer_vector user explicit alpha-reduction with multi-passes
        vector = self._doc2vec.infer_vector(doc_words,
                                            alpha=alpha,
                                            min_alpha=min_alpha,
                                            steps=passes)
        pv.set_vector(vector.tolist())
        pv.save()

        return paper_pk

    def infer_papers(self, paper_pks, **kwargs):
        """Infer model vector for papers
        """
        self.check_active()

        pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                           maxval=len(paper_pks), redirect_stderr=True).start()
        count = 0
        for paper_pk in paper_pks:
            self.infer_paper(paper_pk, **kwargs)
            pbar.update(count)
            count += 1
        # close progress bar
        pbar.finish()

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


class TextField(TimeStampedModel):
    """Store which fields of paper and related data are used in model training
    """

    text_field = models.CharField(max_length=100, choices=FIELDS_FOR_MODEL,
                                  null=False, blank=False)

    def __str__(self):
        return self.text_field


class PaperVectorsManager(models.Manager):

    def create(self, **kwargs):
        # enforce that model is active
        if not kwargs.get('model').is_active:
            raise ValidationError('model {0} is not active'
                             .format(kwargs.get('model').name))

        return super(PaperVectorsManager, self).create(**kwargs)


class PaperVectors(TimeStampedModel):
    """ Table Paper - NLP Model relationship

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


class PaperNeighbors(TimeStampedModel):
    """ Table of papers nearest neighbors"""
    lsh = models.ForeignKey('LSH')

    paper = models.ForeignKey(Paper, related_name='neighbors')

    # Primary keys of the k-nearest neighbors papers
    neighbors = ArrayField(models.IntegerField(null=True),
                           size=settings.NLP_MAX_KNN_NEIGHBORS,
                           null=True, blank=True)

    def __str__(self):
        return '{model_name}/{time_lapse}'.format(
            model_name=self.lsh.model.name,
            time_lapse=self.lsh.time_lapse)

    def set_neighbors(self, vector):
        self.neighbors = pad_neighbors(vector)
        self.save(update_fields=['neighbors'])

    def get_neighbors(self):
        return self.neighbors[:self.lsh.model.size]

    class Meta:
        unique_together = ('lsh', 'paper')


class JournalVectorsManager(models.Manager):

    def create(self, **kwargs):
        # enforce that model is active
        if not kwargs.get('model').is_active:
            raise ValidationError('model {0} is not active'
                             .format(kwargs.get('model').name))

        return super(JournalVectorsManager, self).create(**kwargs)


class JournalVectors(TimeStampedModel):
    """ Table Journal - NLP Model relationship

    Vector field length is fixed to NLP_MAX_VECTOR_SIZE:
    i.e all model must have embedding space < NLP_MAX_VECTOR_SIZE.
    Shorter vectors are pad with None.

    Use set_vector() to pad and set vector list
    """

    journal = models.ForeignKey(Journal, related_name='vectors')

    model = models.ForeignKey(Model)

    vector = ArrayField(models.FloatField(null=True),
                        size=settings.NLP_MAX_VECTOR_SIZE,
                        null=True, db_index=True)

    objects = JournalVectorsManager()

    class Meta:
        unique_together = ('journal', 'model')

    def __str__(self):
        return '{short_title}/{name}'\
            .format(short_title=self.journal.short_title, name=self.model.name)

    def set_vector(self, vector):
        self.vector = pad_vector(vector)
        self.save()

    def get_vector(self):
        return self.vector[:self.model.size]

class LSHManager(models.Manager):

    def create(self, **kwargs):
        """Return LSH instance, updated"""
        obj = super(LSHManager, self).create(**kwargs)
        # Init lsh
        obj.full_update(**kwargs)
        return obj

    def get(self, **kwargs):
        # print('Warning: Not loading LSH object. Use load() instead')
        return super(LSHManager, self).get(**kwargs)

    def load(self, **kwargs):
        """Return LSH instance loaded from file"""
        obj = super(LSHManager, self).get(**kwargs)
        try:
            # if not on volume try download from s3
            if not os.path.isfile(os.path.join(settings.NLP_LSH_PATH,
                                               '{model_name}-tl{time_lapse}.lsh'.format(
                                               model_name=obj.model.name,
                                               time_lapse=obj.time_lapse))):
                obj.download_from_s3()
            obj.lsh = joblib.load(
                os.path.join(settings.NLP_LSH_PATH,
                             '{model_name}-tl{time_lapse}.lsh'.format(
                                 model_name=obj.model.name,
                                 time_lapse=obj.time_lapse)))
        except EnvironmentError:      # OSError or IOError...
            raise

        return obj


class LSH(TimeStampedModel, S3Mixin):
    """Local Sensitive Hashing to retrieve approximate k-neighbors"""

    # For S3 Mixin
    BUCKET_NAME = settings.NLP_LSH_BUCKET_NAME
    PATH = settings.NLP_LSH_PATH
    AWS_ACCESS_KEY_ID = settings.DJANGO_AWS_ACCESS_KEY_ID
    AWS_SECRET_ACCESS_KEY = settings.DJANGO_AWS_SECRET_ACCESS_KEY

    model = models.ForeignKey(Model, related_name='lsh')

    state = models.CharField(default='NON', choices=LSHF_STATES, max_length=3)

    time_lapse = models.IntegerField(default=-1,
                                     choices=NLP_TIME_LAPSE_CHOICES,
                                     verbose_name='Days from right now')

    lsh = MyLSHForest()

    objects = LSHManager()

    class Meta:
        unique_together = ('model', 'time_lapse')

    def __str__(self):
        return '{0}/{1}/{2}'.format(self.id, self.model.name, self.time_lapse)

    def save(self, *args, **kwargs):
        if not self.state == 'BUS':
            if not os.path.exists(settings.NLP_LSH_PATH):
                os.makedirs(settings.NLP_LSH_PATH)
            joblib.dump(self.lsh, os.path.join(settings.NLP_LSH_PATH,
                '{model_name}-tl{time_lapse}.lsh'.format(
                    model_name=self.model.name,
                    time_lapse=self.time_lapse)))
            self.save_db_only()
        else:
            raise InvalidState('LSH state is {0}'.format(self.state))

    def delete(self, *args, **kwargs):
        try:
            os.remove(
                os.path.join(settings.NLP_LSH_PATH,
                             '{model_name}-tl{time_lapse}.lsh'.format(
                                model_name=self.model.name,
                                time_lapse=self.time_lapse)))
        except IOError:
            pass
        super(LSH, self).delete(*args, **kwargs)

    def save_db_only(self, *args, **kwargs):
        super(LSH, self).save(*args, **kwargs)

    def set_state(self, state):
        self.state = state
        self.save_db_only()

    def get_data(self, partial=False):
        """Get and reshape data for training LSH

        Args:
            partial (bool): If true, perform partial fit of LSH (default=False)

        Returns:
            (np.array): 2D array (sample x vector size) of paper vector
            (list): List of PaperVectos primary keys
            (list): List of corresponding Paper primary keys

        """

        logger.info(
            'Updating LSH ({pk}/{model_name}/{time_lapse}) - getting data...'
                .format(pk=self.id, model_name=self.model.name,
                        time_lapse=self.time_lapse))

        # size of the model space used for truncation of vector
        vec_size = self.model.size

        if self.time_lapse > 0:
            from_date = (timezone.now() -
                         timezone.timedelta(
                             days=self.time_lapse)).date()
            data = PaperVectors.objects\
                .filter(Q(model=self.model) &
                        (Q(paper__date_ep__gt=from_date) |
                         (Q(paper__date_pp__gt=from_date) &
                          Q(paper__date_ep=None))))\
                .exclude(Q(paper__is_trusted=False) | Q(paper__abstract=''))\
                .values_list('pk', 'paper__pk', 'vector')
        else:
            if partial:
                data = PaperVectors.objects\
                    .filter(model=self.model, is_in_full_lsh=False)\
                    .exclude(Q(paper__is_trusted=False) | Q(paper__abstract=''))\
                    .values_list('pk', 'paper__pk', 'vector')
            else:
                data = PaperVectors.objects\
                    .filter(model=self.model)\
                    .exclude(Q(paper__is_trusted=False) | Q(paper__abstract=''))\
                    .values_list('pk', 'paper__pk', 'vector')

        # Reshape data
        pv_pks = []
        new_pks = []    # use for partial update
        x_data = np.zeros((data.count(), vec_size))
        for i, dat in enumerate(data):
            # store PaperVector pk
            pv_pks.append(data[i][0])
            # store paper pk
            new_pks.append(data[i][1])
            # build input matrix for fit
            x_data[i, :] = data[i][2][:vec_size]

        # store paper pks in lsh
        self.lsh.pks += new_pks

        return x_data, pv_pks, new_pks

    def update_if_full_lsh_flag(self, pv_pks):
        """Update is_in_full_lsh flag of PaperVectors"""
        # this is a the full LSH
        if self.time_lapse < 0:
            PaperVectors.objects.filter(pk__in=pv_pks)\
                .update(is_in_full_lsh=True)

    def _update(self, partial=False, n_estimators=10, n_candidates=50, **kwargs):
        # TODO: Complete docstring
        """Update LSH with new data and fit

        Partial update is allowed only for LSH with time_lapse=-1. Partial fit
        cannot be performed on LSH with time_lapse > 0 because samples cannot
        be removed from LSH, only added.

        Args:
            partial (bool): if True, do partial update
            n_estimators (int): Number of trees in the LSH Forest
            n_candidates (int): Minimum number of candidates evaluated per
                estimator
        """

        # Register status as busy
        self.set_state('BUS')

        logger.info(
            'Update LSH ({pk}/{model_name}/{time_lapse}) - start...'
                .format(pk=self.id, model_name=self.model.name,
                        time_lapse=self.time_lapse))
        t0 = time.time()

        if not partial:         # Init LSH
            self.lsh = MyLSHForest(n_estimators=n_estimators,
                                   n_candidates=n_candidates)

        # Get data
        data, pv_pks, new_pks = self.get_data(partial=partial)

        # Train
        if pv_pks:      # if data

            logger.info(
                'Update LSH ({pk}/{model_name}/{time_lapse}) - fit...'
                    .format(pk=self.id,
                            model_name=self.model.name,
                            time_lapse=self.time_lapse))

            if partial:
                self.lsh.partial_fit(data)
            else:
                self.lsh.fit(data)

            # Register status as idle
            self.set_state('IDL')

            # Save
            self.save()

            # Update PaperVector.is_in_full_lsh
            self.update_if_full_lsh_flag(pv_pks)

            # Update neighbors
            if partial:
                self.update_neighbors(new_pks)
            else:
                self.update_neighbors()

        tend = time.time() - t0
        logger.info(
            'full update of LSH ({pk}/{model_name}/{time_lapse}) done in {time}'
                .format(pk=self.id, model_name=self.model.name,
                        time_lapse=self.time_lapse, time=tend))

    def update(self, **kwargs):
        """Update dispatcher for LSH

        Depending on time_lapse, update is routed to partial udpate or
        full_update
        """
        # Register state as busy
        self.set_state('BUS')

        if 'partial' in kwargs and kwargs['partial']:
            if self.time_lapse < 0:
                self._update(**kwargs)
            else:
                raise ValueError('partial can be true with lsh time_lapse '
                                 'defined')
        else:
            if self.time_lapse < 0:
                self._update(partial=True, **kwargs)
            else:
                self._update(partial=False, **kwargs)

        # Register state as busy
        self.set_state('IDL')

    def full_update(self, **kwargs):
        """Full update of LSH"""
        # Register state as busy
        self.set_state('BUS')
        self._update(partial=False, **kwargs)
        # Register state as busy
        self.set_state('IDL')

    def update_neighbors(self, *args):
        """Update neighbors of papers

        Args:
            Optional:
            args: list or query set of paper primary keys to update
        """

        # Add Neighbors to PaperNeighbors
        if not args:  # populate all neighbors associated with LSH
            pks = self.lsh.pks
        else:
            pks = args[0]

        for count, pk in enumerate(pks):

            if not count % np.ceil(len(pks)/10):
                logger.info(
                    'Updating LSH ({pk}/{model_name}/{time_lapse}) - updating '
                    'neighbors ({perc:.0f}%)...'.format(
                        pk=self.id,
                        model_name=self.model.name,
                        time_lapse=self.time_lapse,
                        perc=np.round(count / np.ceil(len(pks)/10.) * 10)))

            self.populate_neighbors(pk)
        # self.populate_neighbors_bulk(pks)


    def populate_neighbors(self, paper_pk):
        """Populate neighbors of paper
        """

        pv = PaperVectors.objects.get(paper_id=paper_pk, model=self.model)
        vec = pv.get_vector()

        pks = self.k_neighbors(vec,
                               n_neighbors=settings.NLP_MAX_KNN_NEIGHBORS + 1)

        pks = pks.flatten()[1:]     # remove first element (self)

        pn, _ = PaperNeighbors.objects.get_or_create(lsh_id=self.pk,
                                                     paper_id=paper_pk)
        pn.set_neighbors(pks)


    # BELOW IS A FAILED TENTATIVE TO SPEED UP POPULATE NEIGHBORS
    # def populate_neighbors_bulk(self, paper_pks):
    #     pvs = PaperVectors.objects.filter(paper_id__in=paper_pks, model=self.model)\
    #         .values('pk', 'paper__pk', 'vector')
    #
    #     mat = [pv['vector'][:self.model.size] for pv in pvs]
    #     mat_np = np.array(mat)
    #     indices = self.k_neighbors(mat_np,
    #                                n_neighbors=settings.NLP_MAX_KNN_NEIGHBORS + 1)
    #     ind = {}
    #     for i, pv in enumerate(pvs):
    #         ind[pv['paper__pk']] = indices[i, :][1:].flatten()
    #
    #     pns = PaperNeighbors.objects.filter(lsh_id=self.pk,
    #                                         paper_id__in=paper_pks)
    #     p_pns_pks = [pn.paper.pk for pn in pns]
    #     p_pks_create = [pk for pk in paper_pks if pk not in p_pns_pks]
    #     # bulk create pn
    #     pns_create = []
    #     for pk in p_pks_create:
    #         pns_create.append(PaperNeighbors(lsh_id=self.pk,
    #                                          paper_id=pk))
    #     pns_new = PaperNeighbors.objects.bulk_create(pns_create)
    #     pns = list(pns) + pns_new
    #
    #     # bulk update
    #     for pn in pns:
    #         pn.neighbors = pad_neighbors(ind[pn.paper_id])
    #     bulk_update(pns, update_fields=['neighbors'])


    def k_neighbors(self, seed, **kwargs):
        """Return k nearest neighbors

        Arguments:
            seed (np.array or list): A 1d or 2d array (N samples x model.size)

            Optional:
            n_neighbors (int): number of neighbors

        Returns:
            (np.array): A 2d array of indices of the k-nearest neighbors

        """

        if 'n_neighbors' in kwargs:
            n_neighbors = kwargs['n_neighbors']
        else:
            n_neighbors = settings.NLP_MAX_KNN_NEIGHBORS

        if self.state == 'IDL':
            if isinstance(seed, list):
                seed = np.array(seed).squeeze()
            if seed.ndim == 1:
                assert len(seed) == self.model.size
            else:
                assert seed.shape[1] == self.model.size

            distances, indices = self.lsh.kneighbors(seed,
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

    def k_neighbors_pks(self, seed_pks, model_pk, **kwargs):
        """Return k nearest neighbors

        This method duplicates k_neighbors method because it is used for
        asynchronous task so that only list of PaperVector keys is serialized
        and not the entire 2D array (that by the way crashes celery in my test)

        Arguments:
            seed (list): list of primary PaperVector keys

            Optional:
            n_neighbors (int): number of neighbors

        Returns:
            (list): A nested list of indices of the k-nearest neighbors
        """

        data = PaperVectors.objects\
            .filter(paper__pk__in=seed_pks,
                    model__pk=model_pk)\
            .values('vector')
        model_size = Model.objects.get(pk=model_pk).size

        mat = np.zeros((len(data), model_size), dtype=np.float)
        # populate
        for i, entry in enumerate(data):
            mat[i] = np.array(entry['vector'][:model_size])

        indices = self.k_neighbors(mat, **kwargs)
        return indices.tolist()

    def tasks(self, *args, **kwargs):
        """Use from tasks.py for calling task while object remains in-memory

        This is an odd design. This method dispatches tasks that can be run
        from LSH. It is so because we want to have the lsh object in-memory
        to avoid over-head of loading from file. Given the Task class of Celery
        the work-around I found is to define a BIG task that dispatches to
        smaller tasks. Therefore the instance of the tasks is associated to
        one single thread/worker and everything is 'fast'

        Args:
            task (string): A string defining the task. Either 'update',
                'full_update', 'k_neighbors', 'populate_neighbors'

            Optional Kwargs for:
            'k_neighbors':
                vec (np.array or list): An array corresponding to vector
                n_neighbors (int): An int defining the number of neighbors to
                    search
            'populate_neighbors':
                paper_pk (int): The primary key of a Paper

        Returns:
            Depends on the task
        """

        # the following check is awkward but it is related to
        # https://github.com/celery/celery/issues/2695 and the fact that
        # passing kwargs argument into chain is buggy currently
        if len(args) > 1:  # case: task = 'populate_neighbors'
            task = args[1]
            paper_pk = args[0]
        else:
            task = args[0]

        # update of LSH (partial or not)
        if task == 'update':
            self.update()
            return 0
        # full update of LSH
        elif task == 'full_update':
            self.full_update()
            return 0
        elif task == 'k_neighbors_pks':
            try:
                seed_pks = kwargs['seed_pks']
                model_pk = kwargs['model_pk']
                k = kwargs['k']
            except KeyError as e:
                raise e
            pks = self.k_neighbors_pks(seed_pks, model_pk, n_neighbors=k)
            return pks
        # populate the PaperNeighbors for paper_pk
        elif task == 'populate_neighbors':
            self.populate_neighbors(paper_pk)
            return paper_pk
        else:
            print(task)
            raise ValueError('Unknown task action')