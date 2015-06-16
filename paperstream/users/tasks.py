import logging
from django.conf import settings
from config.celery import celery_app as app

logger = logging.getLogger(__name__)


@app.task(name='tasks.update_lib')
def update_lib(user, provider):

    # get social
    social = user.social_auth.get(provider=provider)

    # get backend
    backend = social.get_backend_instance()

    # build session
    session = backend.get_session(social, user)

    # update lib
    backend.update_lib(user, session)