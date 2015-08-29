# -*- coding: utf-8 -*-

import environ
ROOT_DIR = environ.Path(__file__) - 3
print(type(ROOT_DIR))
print(ROOT_DIR.path('nlp','test.log'))
