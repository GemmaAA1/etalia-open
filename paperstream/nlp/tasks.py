import logging
from config.celery import celery_app as app
from celery import chain, Task
from nlp.models import Model, LSH
from core.constants import NLP_TIME_LAPSE_CHOICES

logger = logging.getLogger(__name__)

class EmbedPaperTask(Task):
    """Abstract Embedding Paper task for model
    Use to load model in __init__ so that it is cached for subsequent call
    to task
    """

    ignore_result = False
    model_name = None
    _model = None

    def __init__(self, *args, **kwargs):
        if 'model_name' in kwargs:
            model_name = kwargs['model_name']
            # check in model_name is known
            choices = [model['name'] for model in
                       Model.objects.all().values('name')]
            if model_name in choices:
                self.model_name = model_name
            else:
                raise ValueError('<model_name> unknown, choices are: {0}'
                                 .format(choices))

    @property
    def model(self):
        if self._model is None:
            self._model = Model.objects.load(name=self.model_name)
            return self._model

        last_modified = Model.objects.get(name=self.model_name).modified
        if not self._model.modified == last_modified:
            self._model = Model.objects.load(name=self.model_name)

        return self._model

    def run(self, paper_pk):
        return self.model.infer_paper(paper_pk)


class LSHTask(Task):
    """Abstract LSH task for model
    Use to load lsh instance in __init__ so that it is cached for subsequent call
    to task
    """

    ignore_result = False
    model_name = None
    time_lapse = None
    _lsh = None

    def __init__(self, *args, **kwargs):
        if ('model_name' in kwargs) and ('time_lapse' in kwargs):
            model_name = kwargs['model_name']
            time_lapse = kwargs['time_lapse']

            # check if model_name is known
            choices = [model['name'] for model in
                       Model.objects.all().values('name')]
            if model_name in choices:
                self.model_name = model_name
            else:
                raise ValueError('<model_name> unknown, choices are: {0}'
                                 .format(choices))
            # check if time_lapse allowed
            choices = [time_lapse for time_lapse, _ in NLP_TIME_LAPSE_CHOICES] + \
                      [None]
            if time_lapse in choices:
                self.time_lapse = time_lapse
            else:
                raise ValueError('<{0}> not in possible time_lapse choices'
                                 .format(time_lapse))

    @property
    def lsh(self):
        if self._lsh is None:
            self._lsh = LSH.objects.load(model__name=self.model_name,
                                         time_lapse=self.time_lapse)
            return self._lsh
        # if lsh has been modified, reload
        last_modified = LSH.objects.get(model__name=self.model_name,
                                        time_lapse=self.time_lapse).modified
        if not self._lsh.modified == last_modified:
            self._lsh = LSH.objects.load(model__name=self.model_name,
                                         time_lapse=self.time_lapse)
        return self._lsh

    def run(self, *args, **kwargs):
        return self.lsh.tasks(*args, **kwargs)

# Model based tasks factory
def register_all_models_and_lshs_tasks():
    # Create embedding task from model
    model_names = Model.objects.all().values_list('name', flat=True)

    for model_name in model_names:
        cls = EmbedPaperTask(model_name=model_name)
        app.task(cls, name='nlp.tasks.embed_paper_{model_name}'.format(
            model_name=model_name))

    # Create lsh related tasks
    for model_name in model_names:
        for time_lapse, _ in NLP_TIME_LAPSE_CHOICES:
            cls = LSHTask(model_name=model_name, time_lapse=time_lapse,
                          bind=True)
            app.task(cls, name='nlp.tasks.lsh_{model_name}_{time_lapse}'
                     .format(model_name=model_name,
                             time_lapse=time_lapse))

register_all_models_and_lshs_tasks()


def embed_all_models(paper_pk):
    """Apply all possible embedding from model for paper_pk
    """
    model_names = Model.objects.all().values_list('name', flat=True)
    for model_name in model_names:
        try:
            embed_task = app.tasks['nlp.tasks.embed_paper_{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Embeding task for {model_name} not defined'.format(
                model_name=model_name))
            continue
        embed_task.apply_async(args=(paper_pk,))


def embed_all_models_and_find_neighbors(paper_pk):
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
                lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_full'.format(
                    model_name=model_name,
                    time_lapse=time_lapse)]
            except KeyError:
                logger.error('LSH task for {model_name}/{time_lapse} not defined'
                             .format(model_name=model_name,
                                     time_lapse=time_lapse))
                continue

            chain(embed_task.s(paper_pk),
                  lsh_task.s(kwargs={'task': 'populate_neighbors'}))