from django.shortcuts import render, redirect
from django.views.generic import UpdateView
from .forms import SocialPrimaryForm
# Create your views here.

def done(request):
    redirect('landing')

def require_primary(request):
    backend = request.session['partial_pipeline']['backend']
    details = request.session['partial_pipeline']['kwargs']['details']
    context = {
        'primary': True,
        'backend': backend,
        'first_name': details.get('first_name', ''),
        'last_name': details.get('last_name', ''),
        'email': details.get('email', ''),
    }
    return render(request, 'users/info.html', context)

# TODO change this into Form generic view
def require_affiliation(request):
    backend = request.session['partial_pipeline']['backend']
    details = request.session['partial_pipeline']['kwargs']['details']
    affiliation = details.get('affiliation', {})
    context = {
        'affiliation': True,
        'backend': backend,
        'department': affiliation.get('department', ''),
        'institution': affiliation.get('institution', ''),
        'city': affiliation.get('city', ''),
        'state': affiliation.get('state', ''),
        'country': affiliation.get('country', ''),
    }
    return render(request, 'users/info.html', context)

def validation_sent(request):
    context = {
        'validation_sent': True,
        'email': request.session.get('email_validation_address')
    }
    return render(request, 'users/info.html', context)

