from django.shortcuts import redirect, render
from django.utils import timezone
from social.pipeline.partial import partial
from .models import Affiliation


@partial
def require_primary(strategy, details, user=None, is_new=False, *args, **kwargs):
    """ Redirect to primary info form for user to check them
    """
    if user and user.email and user.first_name and user.last_name:
        return
    elif is_new:
        email = strategy.request_data().get('email')
        first_name = strategy.request_data().get('first_name')
        last_name = strategy.request_data().get('last_name')
        if email and first_name and last_name:
            details['email'] = email
            details['first_name'] = first_name
            details['last_name'] = last_name
        else:
            return redirect('users:require_primary')

@partial
def require_affiliation(strategy, details, request=None, user=None, *args, **kwargs):
    if getattr(user, 'affiliation'):
        return
    elif request.get('city', None) and request.get('country', None) and user:
        affiliation, _ = Affiliation.objects.get_or_create(
            department=request.get('department', ''),
            institution=request.get('institution', ''),
            city=request.get('city', ''),
            state=request.get('state', ''),
            country=request.get('country', ''),
        )
        user.affiliation = affiliation
        user.save()
        return
    else:
        return redirect('users:require_affiliation')

@partial
def update_user_lib(backend, social, user, *args, **kwargs):
    session = backend.get_session(social, user, *args, **kwargs)
    backend.update_lib(session, user, *args, **kwargs)

    return {}


