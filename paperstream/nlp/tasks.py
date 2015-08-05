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
        # super(EmbedPaperTask, self).__init__(*args, **kwargs)
        if 'model_name' in kwargs:
            self.model_name=kwargs['model_name']

    @property
    def model(self):
        if self._model is None:
            # self._model = Model.objects.load(name='dbow')
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
        if 'model_name' in kwargs:
            self.model_name = kwargs['model_name']
        if 'time_lapse' in kwargs:
            self.time_lapse = kwargs['time_lapse']

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
        cls = LSHTask(model_name=model_name, bind=True)
        app.task(cls, name='nlp.tasks.lsh_{model_name}_full'
                 .format(model_name=model_name, time_lapse=None))
        # time_lapse dependant LSHs
        for time_lapse in TIME_LAPSE_CHOICES:
            cls = LSHTask(model_name=model_name, bind=True)
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
        # Send task for LSH full
        try:
            embed_task = app.tasks['nlp.tasks.embed_paper_{model_name}'.format(
                model_name=model_name)]
        except KeyError:
            logger.error('Embeding task for {model_name} not defined'.format(
                model_name=model_name))
            continue
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
                    time_lapse=time_lapse)]
            except KeyError:
                logger.error('LSH task for {model_name}/{time_lapse} not defined'
                             .format(model_name=model_name,
                                     time_lapse=time_lapse))
                continue

            chain(embed_task.s(),
                  lsh_task.s(kwargs={'task': 'populate_neighbors'}))
            chain.delay(paper_pk)
