import os
from django.conf import settings
from django.apps import AppConfig

class NLPConfig(AppConfig):
    """Use to create default folders that stores nlp models and data
    """

    name = 'paperstream.nlp'

    # Creating default folders for NLP
    if not(os.path.isdir(settings.NLP_DOC2VEC_PATH)):
        os.makedirs(settings.NLP_DOC2VEC_PATH, exist_ok=True)
    if not(os.path.isdir(settings.NLP_DATA_PATH)):
        os.makedirs(settings.NLP_DATA_PATH, exist_ok=True)
    if not(os.path.isdir(settings.NLP_LSH_PATH)):
        os.makedirs(settings.NLP_LSH_PATH, exist_ok=True)
