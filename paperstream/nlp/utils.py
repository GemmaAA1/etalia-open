# To be run from paperstream shell
import os
import re
import nltk
import glob
from bs4 import BeautifulSoup

from gensim.models import Doc2Vec
from gensim.models import Phrases
from gensim.models.doc2vec import TaggedDocument
from sklearn.neighbors import LSHForest

import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)

def paper2tokens(paper, **kwargs):
    """Convert paper instance to tokens list
    """
    if 'fields' in kwargs:
        if not isinstance(kwargs['fields'], list):
            raise TypeError('<fields> must be list of field strings')
        fields = kwargs['fields']
    else:
        raise ValueError('<fields> must be in kwargs')

    tokens_list = []
    for j, field in enumerate(fields):
        # add trailling '. ' after field if no punctuation mark found
        field_val = getattr(paper, field)
        if field_val:
            if not field_val.strip()[-1] in ['.', '?', '!']:
                str_val = '{0}. '.format(field_val.strip())
            else:
                str_val = '{0} '.format(field_val.strip())
        else:
            str_val = ''

        # pre-process
        tokens_list.append(pre_process_text(str_val))

    return ' '.join(tokens_list).split(' ')

def pre_process_text(text):
    """Preprocess text
    """
    # remove HTML
    text = BeautifulSoup(text).get_text()

    # lower text
    text = text.lower()

    # tokenize
    tokens = nltk.word_tokenize(text)

    return ' '.join(tokens)

class TaggedDocumentsIterator(object):
    """Iterator of TaggedDocument
    """
    def __init__(self, dir_path, **kwargs):
        self.FILE_FORMAT = '*.txt'
        if 'phraser' in kwargs:
            if not isinstance(kwargs['phraser'], Phrases):
                raise TypeError('phraser not a Phrases instance')
            else:

                self.phraser = kwargs['phraser']
        else:
            self.phraser = None
        self.dir_path = dir_path
        self.filenames = glob.glob(os.path.join(dir_path, self.FILE_FORMAT))

    def __iter__(self):
        for filename in self.filenames:
            for line in open(filename):
                # print(line)
                pk, j_pk, text = re.match(r'(?P<pk>[\d]+), (?P<j_pk>j_[\d]+): (.+)', line).groups()
                if self.phraser:
                    text_l = self.phraser[text.strip().split(' ')]
                else:
                    text_l = text.strip().split(' ')
                yield TaggedDocument(text_l, [pk, j_pk])


class MyLSHForest(LSHForest):
    """Adding pks attribute to store Paper pk correspondence in LSH
    """

    def __init__(self, *args, **kwargs):
        super(MyLSHForest, self).__init__(*args, **kwargs)
        self.pks = []


def model_attr_getter(attr):
    def get_any(self):
        return getattr(self, attr)
    return get_any

def model_attr_setter(attr, attr2):
    def set_any(self, value):
        setattr(self, attr, value)
        if attr2 == 'dm':
            setattr(self._doc2vec, 'sg',  (value + 1) % 2)
        else:
            setattr(self._doc2vec, attr2, value)
    return set_any


class MyDoc2Vec(Doc2Vec):

    def infer_vector_with_seed(self, doc_words, alpha=0.1, min_alpha=0.0001,
                               steps=5, seed=None):
        """Infer a vector for given post-bulk training document.

        Document should be a list of (word) tokens.
        And seed an array of size vector_size
        """

        doctag_vectors = empty((1, self.vector_size), dtype=REAL)
        # doctag_vectors[0] = self.seeded_vector(' '.join(doc_words))
        doctag_vectors[0] = seed
        doctag_locks = ones(1, dtype=REAL)
        doctag_indexes = [0]

        work = zeros(self.layer1_size, dtype=REAL)
        if not self.sg:
            neu1 = matutils.zeros_aligned(self.layer1_size, dtype=REAL)

        for i in range(steps):
            if self.sg:
                train_document_dbow(self, doc_words, doctag_indexes, alpha, work,
                                    learn_words=False, learn_hidden=False,
                                    doctag_vectors=doctag_vectors, doctag_locks=doctag_locks)
            elif self.dm_concat:
                train_document_dm_concat(self, doc_words, doctag_indexes, alpha, work, neu1,
                                         learn_words=False, learn_hidden=False,
                                         doctag_vectors=doctag_vectors, doctag_locks=doctag_locks)
            else:
                train_document_dm(self, doc_words, doctag_indexes, alpha, work, neu1,
                                  learn_words=False, learn_hidden=False,
                                  doctag_vectors=doctag_vectors, doctag_locks=doctag_locks)
            alpha = ((alpha - min_alpha) / (steps - i)) + min_alpha

        return doctag_vectors[0]