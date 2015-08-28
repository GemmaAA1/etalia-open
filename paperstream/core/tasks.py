import logging
from celery.canvas import chain
from config.celery import celery_app as app
from paperstream.core.constants import NLP_TIME_LAPSE_CHOICES
from paperstream.nlp.models import Model

logger = logging.getLogger(__name__)


def embed_all_models_and_find_neighbors(paper_pk):
    """
    """
    model_names = Model.objects.all().values_list('name', flat=True)
    for model_name in model_names:
        # Send task for embedding
        try:
            embed_task = app.tasks['nlp.tasks.embed_paper_{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Embeding task for {model_name} not defined'.format(
                model_name=model_name))
            continue

        # Send task for time_lapse related LSHs
        for time_lapse, _ in NLP_TIME_LAPSE_CHOICES:
            try:
                lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_{time_lapse}'
                    .format(model_name=model_name, time_lapse=time_lapse)]
            except KeyError:
                logger.error('LSH task for {model_name}/{time_lapse} not defined'
                             .format(model_name=model_name,
                                     time_lapse=time_lapse))
                continue

            task = chain(embed_task.s(paper_pk),
                         lsh_task.s('populate_neighbors'))
            task.delay()
