import os
from django.conf import settings
from django.apps import AppConfig

from paperstream.core.utils import makedirs_p

class NLPConfig(AppConfig):
    """Use to create default folders that stores nlp models and data
    """

    name = 'paperstream.nlp'

    # Creating default folders for NLP
    if not(os.path.isdir(settings.NLP_DOC2VEC_PATH)):
        makedirs_p(settings.NLP_DOC2VEC_PATH)
    if not(os.path.isdir(settings.NLP_DATA_PATH)):
        makedirs_p(settings.NLP_DATA_PATH)
    if not(os.path.isdir(settings.NLP_LSH_PATH)):
        makedirs_p(settings.NLP_LSH_PATH)



