from config.celery import celery_app as app
from .models import ConsumerPubmed, ConsumerArxiv, ConsumerElsevier

@app.task(name='tasks.pubmed_run')
def pubmed_run(self, name):
    pubmed_consumer = ConsumerPubmed.object.get(name=name)
    pubmed_consumer.run_once_per_period()

@app.task(name='tasks.arxiv_run')
def arxiv_run(self, name):
    arxiv_consumer = ConsumerArxiv.object.get(name=name)
    arxiv_consumer.run_once_per_period()

@app.task(name='tasks.elsevier_run')
def elsevier_run(self, name):
    elsevier_consumer = ConsumerElsevier.object.get(name=name)
    elsevier_consumer.run_once_per_period()
