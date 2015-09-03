# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import, print_function

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

from django.conf import settings

import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)


class MyDoc2Vec(Doc2Vec):
    pass

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
    # add data to path
    nltk.data.path += [settings.NLP_NLTK_DATA_PATH]

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
            with open(filename) as fp:
                for line in fp:
                    pk, j_pk, text = re.match(
                        r'(?P<pk>[\d]+), (?P<j_pk>j_[\d]+): (.+)', line)\
                        .groups()
                    if self.phraser:
                        text_l = self.phraser[text.strip().split(u' ')]
                    else:
                        text_l = text.strip().split(u' ')
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