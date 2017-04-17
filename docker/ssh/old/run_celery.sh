#!/bin/sh

# Wait for worker to be up
wait 5

# run Celery worker for our project myproject with Celery configuration stored in Celeryconf
su -m myuser -c "celery worker -A config -Q default,nlp,pe,te,udpate_engines,feed,consumers,altmetric,library,beat,broadcast_pe_tasks,broadcast_te_tasks,test --config config.celery_settings.development -n default@%h --loglevel=INFO"