import logging
from config.celery import celery_app as app
from users.backends.mendeley import CustomMendeleyOAuth2

logger = logging.getLogger(__name__)
