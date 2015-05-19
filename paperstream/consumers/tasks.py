from config.celery import celery_app as app
from .models import ConsumerPubmed

@app.task(name='tasks.pubmed_once_a_day')
def pubmed_run_once_a_day(self, name):
    pubmed_consumer = ConsumerPubmed.object.get(name=name)
    pubmed_consumer.run_once_a_day()
