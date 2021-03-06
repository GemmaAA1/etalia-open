# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
from django.conf import settings
from django.apps import AppConfig

from etalia.core.utils import makedirs_p


class NLPConfig(AppConfig):
    """Use to create default folders that stores nlp models and data
    """

    name = 'etalia.nlp'
    verbose_name = "NLP app"

    # Creating default folders for NLP
    if not(os.path.isdir(settings.NLP_MODELS_PATH)):
        makedirs_p(settings.NLP_MODELS_PATH)
    if not(os.path.isdir(settings.NLP_DATA_PATH)):
        makedirs_p(settings.NLP_DATA_PATH)
    if not(os.path.isdir(settings.NLP_PE_PATH)):
        makedirs_p(settings.NLP_PE_PATH)
    if not(os.path.isdir(settings.NLP_TE_PATH)):
        makedirs_p(settings.NLP_TE_PATH)

    def ready(self):
        import etalia.nlp.signals  #noqa
