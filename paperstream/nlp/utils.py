# To be run from paperstream shell
import os
import re
import nltk
import glob
from bs4 import BeautifulSoup
from django.conf import settings

from gensim.models.doc2vec import TaggedDocument

import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)


class TaggedDocumentsIterator(object):

    def __init__(self, dir_path):
        self.FILE_FORMAT = '*.txt'
        self.dir_path = dir_path
        self.filenames = glob.glob(os.path.join(dir_path, self.FILE_FORMAT))

    def __iter__(self):
        #TODO: shuffle lines ?
        for filename in self.filenames:
            for line in open(filename):
                # print(line)
                pk, text = re.match(r'([\d]+): (.+)', line).groups()
                yield TaggedDocument(text.strip().split(' '), [pk, ])


class DumpPaperData(object):
    """Dump papers data to pre-process text files.
    Papers data are separated in multiple file to spare memory when building
    model. File are composed of N (CHUNK_SIZE) documents. Each line of the
    files start with the primary key of the document and then the pre-process
    string corresponding to the document field(s) dumped.
     ex:
        123: the title of the document #1 . the abstract of the document .
        124: the title of the document #2 . the abstract of the document .
    """

    def __init__(self, **kwargs):

        self.CHUNK_SIZE = settings.NLP_CHUNK_SIZE

        if 'fields' in kwargs:
            if not isinstance(kwargs['fields'], list):
                raise TypeError('<fields> must be list of field strings')
        if 'to' in kwargs:
            if not isinstance(kwargs['to'], str):
                raise TypeError('<to> must be path string')


        self.fields = kwargs.get('fields', ['abstract'])
        self.to = kwargs.get('to', settings.NLP_DATA_PATH)

    @staticmethod
    def pre_process_text(text):
        # remove HTML
        text = BeautifulSoup(text).get_text()

        # lower text
        text = text.lower()

        # tokenize
        tokens = nltk.word_tokenize(text)

        return ' '.join(tokens)

    def dump(self, papers):

        tot = papers.count()
        file_count = 0
        file = None

        if not os.path.exists(self.to):
            os.makedirs(self.to)

        for count, paper in enumerate(papers):

            if not count % 10:
                print('{0}/{1}'.format(count, tot))

            if not count % self.CHUNK_SIZE:
                if file:
                    file.close()
                file = open(os.path.join(self.to,
                                         '{0:06d}.txt'.format(file_count)),
                            'w+')
                file_count += 1

            # write head
            file.write('{pk}: '.format(pk=paper.pk))

            # write line body
            line_val = []
            for j, field in enumerate(self.fields):
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
                str_val = self.pre_process_text(str_val)

                line_val.append(str_val)

            # write to file
            file.write(' '.join(line_val).strip())

            # write new line
            file.write('\n')

        file.close()


