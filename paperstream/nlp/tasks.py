import logging
from config.celery import celery_app as app
from celery import chain
from nlp.models import Model, LSH
from core.constants import TIME_LAPSE_CHOICES

logger = logging.getLogger(__name__)

class EmbedPaperTask(app.Task):
    """Abstract Embedding Paper task for model
    Use to load model in __init__ so that it is cached for subsequent call
    to task
    """
    abstract = True
    ignore_result = True
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
        else:
            raise TypeError('missing <model_name> argument')

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


class LSHTask(app.Task):
    """Abstract LSH task for model
    Use to load lsh instance in __init__ so that it is cached for subsequent call
    to task
    """
    abstract = True
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
            choices = [time_lapse[0] for time_lapse in TIME_LAPSE_CHOICES] + \
                      [None]
            if time_lapse in choices:
                self.time_lapse = time_lapse
            else:
                raise ValueError('<{0}> not in possible time_lapse choices'
                                 .format(time_lapse))
        else:
            raise TypeError('<model_name> and <time_lapse> must be defined')

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
def register_model_tasks():
    # Create embedding task from model
    model_names = Model.objects.all().values_list('name', flat=True)

    for model_name in model_names:
        cls = EmbedPaperTask(model_name=model_name, bind=True)
        app.task(cls, name='nlp.tasks.embed_paper_{model_name}'.format(
            model_name=model_name))

    # Create lsh related tasks
    for model_name in model_names:
        # Full LSH
        cls = LSHTask(model_name=model_name, time_lapse=None, bind=True)
        app.task(cls, name='nlp.tasks.lsh_{model_name}_full'
                 .format(model_name=model_name, time_lapse=None))
        # time_lapse dependant LSHs
        for time_lapse in TIME_LAPSE_CHOICES:
            cls = LSHTask(model_name=model_name, time_lapse=time_lapse[0],
                          bind=True)
            app.task(cls, name='nlp.tasks.lsh_{model_name}_{time_lapse}'
                     .format(model_name=model_name,
                             time_lapse=time_lapse[0]))
register_model_tasks()


def all_embedings_paper(paper_pk):
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


def all_embeddings_and_neighbors(paper_pk):
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

        # Send task for LSH full
        try:
            lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_full'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('LSH task for {model_name}/full not defined'.format(
                model_name=model_name))
            continue

        chain(embed_task.s(),
              lsh_task.s(kwargs={'task': 'populate_neighbors'}))
        chain.delay(paper_pk)

        # Send task for time_lapse related LSHs
        for time_lapse in TIME_LAPSE_CHOICES:
            try:
                lsh_task = app.tasks['nlp.tasks.lsh_{model_name}_full'.format(
                    model_name=model_name,
                    time_lapse=time_lapse[0])]
            except KeyError:
                logger.error('LSH task for {model_name}/{time_lapse} not defined'
                             .format(model_name=model_name,
                                     time_lapse=time_lapse[0]))
                continue

            chain(embed_task.s(),
                  lsh_task.s(kwargs={'task': 'populate_neighbors'}))
            chain.delay(paper_pk)
