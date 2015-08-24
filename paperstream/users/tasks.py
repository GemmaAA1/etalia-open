import logging
from django.contrib.auth import get_user_model
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

User = get_user_model()


@app.task()
def update_lib(user_pk, provider_name):
    """Async task for updating user library"""
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

    return user_pk
