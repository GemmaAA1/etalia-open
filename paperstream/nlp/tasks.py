from config.celery import celery_app as app
from nlp.models import Model, LSH

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


class LSHTask(app.Task):
    """Abstract LSH task for model
    Use to load lsh instance in __init__ so that it is cached for subsequent call
    to task
    """
    abstract = True
    ignore_result = False
    model_name = None
    _lsh = None

    def __init__(self, *args, **kwargs):
        if 'model_name' in kwargs:
            self.model_name=kwargs['model_name']

    @property
    def lsh(self):
        if self._lsh is None:
            self._lsh = LSH.objects.load(model__name=self.model_name)
            return self._lsh
        # if lsh has been modified, reload
        last_modified = LSH.objects.get(model__name=self.model_name).modified
        if not self._lsh.modified == last_modified:
            self._lsh = LSH.objects.load(model__name=self.model_name)
        return self._lsh

# Create embedding task here:
# embed paper for model dbow
@app.task(base=EmbedPaperTask, model_name='dbow', bind=True)
def dbow_embed_paper(paper_pk):
    dbow_embed_paper.model.infer_paper(paper_pk)

# Create lsh related task here:
@app.task(base=LSHTask, model_name='dbow', bind=True)
def dbow_lsh(*args, **kwargs):
    return dbow_lsh.lsh.mono_task(*args, **kwargs)


