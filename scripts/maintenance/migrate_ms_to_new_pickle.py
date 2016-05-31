#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from sklearn.externals import joblib
import os
import sys
import pickle


def setup_django():
    # search for config/ is parent tree directory and setup django
    cpath = os.path.abspath(__file__)
    while cpath and not os.path.isdir('config'):
        os.chdir('..')
    sys.path.append(os.path.abspath(os.getcwd()))
    import django
    django.setup()


def load_old(obj):
    try:
        # if not on volume try download from s3
        if not os.path.isfile(os.path.join(settings.NLP_MS_PATH,
                                           '{0}.ms_data'.format(obj.name))):
            if getattr(settings, 'NLP_MS_BUCKET_NAME', ''):
                obj.pull_from_s3()

        obj.data = joblib.load(os.path.join(settings.NLP_MS_PATH,
                                            '{0}.ms_data'.format(obj.name)))
        obj.index2pk = joblib.load(os.path.join(settings.NLP_MS_PATH,
                                            '{0}.ms_ind2pk'.format(obj.name)))
        obj.index2journalpk = joblib.load(os.path.join(settings.NLP_MS_PATH,
                                            '{0}.ms_ind2journalpk'.format(obj.name)))
        obj.date = joblib.load(os.path.join(settings.NLP_MS_PATH,
                                            '{0}.ms_date'.format(obj.name)))
    except EnvironmentError:      # OSError or IOError...
        raise

    return obj


def save_new(obj):
    # save files to local volume
    if not os.path.exists(settings.NLP_MS_PATH):
        os.makedirs(settings.NLP_MS_PATH)
    # use joblib for numpy array pickling optimization
    joblib.dump(obj.data, os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_data'))
    with open(os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_ind2pk'), 'wb') as f:
        pickle.dump(obj.index2pk, f)
    with open(os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_ind2journalpk'), 'wb') as f:
        pickle.dump(obj.index2journalpk, f)
    with open(os.path.join(settings.NLP_MS_PATH, obj.name + '.ms_date'), 'wb') as f:
        pickle.dump(obj.date, f)

    # push files to s3
    if obj.BUCKET_NAME:
        obj.upload_state = 'ING'
        obj.save_db_only()
        obj.push_to_s3(ext='pe')
        obj.upload_state = 'IDL'


if __name__ == '__main__':

    setup_django()

    from django.conf import settings
    from etalia.nlp.models import PaperEngine

    mss = PaperEngine.objects.all()
    for ms in mss:
        ms = load_old(ms)
        save_new(ms)


