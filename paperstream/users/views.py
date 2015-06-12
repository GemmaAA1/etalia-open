from django.shortcuts import render, redirect
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView
from django.views.generic.list import ListView
from django.views.generic.detail import DetailView
from django.conf import settings
from django.contrib.auth import get_user_model
from .forms import SocialPrimaryForm
from library.models import Paper, Journal

User = get_user_model()

def logout(request):
    """Logs out user"""
    auth_logout(request)
    return redirect('/')

def done(request):
    redirect('landing')

def require_primary(request):
    backend = request.session['partial_pipeline']['backend']
    details = request.session['partial_pipeline']['kwargs']['details']
    context = {
        'stage': 'primary',
        'backend': backend,
        'first_name': details.get('first_name', ''),
        'last_name': details.get('last_name', ''),
        'email': details.get('email', ''),
    }
    return render(request, 'user/info.html', context)

def require_affiliation(request):
    backend = request.session['partial_pipeline']['backend']
    details = request.session['partial_pipeline']['kwargs']['details']
    affiliation = details.get('tmp_affiliation', {})
    context = {
        'stage': 'affiliation',
        'backend': backend,
        'department': affiliation.get('department', ''),
        'institution': affiliation.get('institution', ''),
        'city': affiliation.get('city', ''),
        'state': affiliation.get('state', ''),
        'country': affiliation.get('country', ''),
    }
    return render(request, 'user/info.html', context)

def validation_sent(request):
    context = {
        'stage': 'validation_sent',
        'email': request.session.get('email_validation_address')
    }
    return render(request, 'user/info.html', context)


class PaperListView(ListView):
    model = Paper
    template_name = 'user/library.html'
    paginate_by = settings.ITEMS_PER_PAGE
    context_object_name = 'paper_list'

    def get_context_data(self, **kwargs):
        context = super(PaperListView, self).get_context_data(**kwargs)
        return context

    def get_queryset(self):
        papers = self.request.user.lib.papers.all()
        return papers

library = PaperListView.as_view()

class ProfileView(DetailView):
    model = User
    template_name = 'user/profile.html'
    fields = ['first_name', 'last_name', 'affiliation']

    def get_object(self, queryset=None):
        return self.request.user

profile = ProfileView.as_view()
