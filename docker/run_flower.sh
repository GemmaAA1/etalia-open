#!/bin/sh

# Wait for worker to be up
wait 10

# run Celery flower for task monitoring
su -m myuser -c "celery flower -A config --config config.celery_settings.development --basic_auth=admin:admin"