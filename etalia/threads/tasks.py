# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

# from config.celery import celery_app as app
#
# logger = logging.getLogger(__name__)
#
# def embed_threads(pks, model_name, batch_size=1000):
#     try:
#         model_task = app.tasks['etalia.nlp.tasks.{model_name}'.format(
#             model_name=model_name)]
#     except KeyError:
#         logger.error('Embeding task for {model_name} not defined'.format(
#             model_name=model_name))
#         raise KeyError
#     pks = list(pks)
#     nb_papers = len(pks)
#     nb_batches = nb_papers // batch_size
#     pks_batched = [pks[i*batch_size:(1+i)*batch_size] for i in range(nb_batches)]
#     pks_batched.append(pks[nb_batches * batch_size:])
#
#     for batch in pks_batched:
#         model_task.delay('infer_threads', batch)