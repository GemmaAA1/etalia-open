#!/bin/sh

# Wait for migrations to be done
wait 10

# run Celery worker for our project myproject with Celery configuration stored in Celeryconf
su -m myuser -c "celery worker -A config -Q default,nlp,pe,te,udpate_engines,feed,consumers,altmetric,library,beat,test --config config.celery_settings.development -n default@%h --loglevel=info"