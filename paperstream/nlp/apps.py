import os
from django.conf import settings
from django.apps import AppConfig

class NLPConfig(AppConfig):
    """Use to create default folders that stores nlp models and data
    """

    name = 'nlp'
    
    # Creating default folders for NLP
    if not(os.path.isdir(settings.NLP_DOC2VEC_PATH)):
        os.mkdir(settings.NLP_DOC2VEC_PATH)
    if not(os.path.isdir(settings.NLP_DATA_PATH)):
        os.mkdir(settings.NLP_DATA_PATH)
    if not(os.path.isdir(settings.NLP_LSH_PATH)):
        os.mkdir(settings.NLP_LSH_PATH)
