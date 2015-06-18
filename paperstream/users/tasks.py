import logging
from django.conf import settings
from config.celery import celery_app as app
from django.contrib.auth import get_user_model
logger = logging.getLogger(__name__)

User = get_user_model()

@app.task(name='tasks.update_lib')
def update_lib(user_pk, provider_name):

    # get user
    user = User.objects.get(pk=user_pk)

    # get social
    social = user.social_auth.get(provider=provider_name)

    # get backend
    backend = social.get_backend_instance()

    # build session
    session = backend.get_session(social, user)

    # update lib
    backend.update_lib(user, session)