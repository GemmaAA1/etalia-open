# To be run from paperstream shell
import os
import re
import nltk
import glob
from progressbar import ProgressBar, Percentage, Bar, ETA
from bs4 import BeautifulSoup

from gensim.models import Phrases
from gensim.models.doc2vec import TaggedDocument

import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)

def paper2tokens(paper, **kwargs):
    """Convert paper instance to tokens list
    """
    if 'fields' in kwargs:
        if not isinstance(kwargs['fields'], list):
            raise TypeError('<fields> must be list of field strings')
    fields = kwargs.get('fields', ['title', 'abstract'])

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

class WordListIterator(object):

    def __init__(self, dir_path):
        self.FILE_FORMAT = '*.txt'
        self.dir_path = dir_path
        self.filenames = glob.glob(os.path.join(dir_path, self.FILE_FORMAT))

    def __iter__(self):
        for filename in self.filenames:
            for line in open(filename):
                # print(line)
                pk, text = re.match(r'([\d]+): (.+)', line).groups()
                yield text.strip().split(' ')


