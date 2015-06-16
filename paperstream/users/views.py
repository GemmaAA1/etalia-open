import json
from django.contrib import messages
from django.shortcuts import HttpResponse
from django.core.urlresolvers import reverse
from django.http import JsonResponse
from django.shortcuts import render, redirect
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView, FormView
from django.views.generic.list import ListView
from django.conf import settings
from django.contrib.auth import get_user_model
from .forms import UserBasicForm, UserAffiliationForm, UserBasicNoEmailForm
from .models import Affiliation
from library.models import Paper

from braces.views import LoginRequiredMixin

User = get_user_model()

def logout(request):
    """Logs out user"""
    auth_logout(request)
    return redirect('/')


def done(request):
    redirect('landing')


class UserBasicInfoSignupView(FormView):

    form_class = UserBasicForm
    template_name = 'user/info.html'

    def get_success_url(self, **kwargs):
        return reverse('social:complete',
        kwargs={'backend': self.request.session['partial_pipeline']['backend']})

    def get_context_data(self, **kwargs):
        context = super(UserBasicInfoSignupView, self).get_context_data(**kwargs)
        context['stage'] = 'basic_info'
        return context

    def get_initial(self):
        initial = super(UserBasicInfoSignupView, self).get_initial()
        details = self.request.session['partial_pipeline']['kwargs']['details']
        initial['first_name'] = details['first_name']
        initial['last_name'] = details['last_name']
        initial['email'] = details['email']
        return initial

    def form_valid(self, form):
        self.request.session['basic_info'] = {
            'first_name': form.cleaned_data['first_name'],
            'last_name': form.cleaned_data['last_name'],
            'email': form.cleaned_data['email']}
        return super(UserBasicInfoSignupView, self).form_valid(form)

require_basic_info = UserBasicInfoSignupView.as_view()


class UserAffiliationSignupView(FormView):

    form_class = UserAffiliationForm
    template_name = 'user/info.html'

    def get_success_url(self, **kwargs):
        return reverse('social:complete',
        kwargs={'backend': self.request.session['partial_pipeline']['backend']})

    def get_context_data(self, **kwargs):
        context = super(UserAffiliationSignupView, self).get_context_data(**kwargs)
        context['stage'] = 'affiliation'
        return context

    def get_initial(self):
        initial = super(UserAffiliationSignupView, self).get_initial()
        details = self.request.session['partial_pipeline']['kwargs']['details']
        return dict(details['tmp_affiliation'], **initial)

    def form_invalid(self, form):
        if Affiliation.objects.filter(**form.cleaned_data).exists():  # invalid because of uniqueness
            affiliation = Affiliation.objects.get(**form.cleaned_data)
            self.request.session['affiliation_pk'] = affiliation.pk
            return super(UserAffiliationSignupView, self).form_valid(form)
        else:
            return super(UserAffiliationSignupView, self).form_invalid(form)

    def form_valid(self, form):
        affiliation = form.save()
        self.request.session['affiliation_pk'] = affiliation.pk
        return super(UserAffiliationSignupView, self).form_valid(form)

require_affiliation = UserAffiliationSignupView.as_view()


def validation_sent(request):
    context = {
        'stage': 'validation_sent',
        'email': request.session.get('email_validation_address')
    }
    return render(request, 'user/info.html', context)


class UserLibraryView(LoginRequiredMixin, ListView):
    model = Paper
    template_name = 'user/library.html'
    paginate_by = settings.ITEMS_PER_PAGE
    context_object_name = 'paper_list'

    def get_context_data(self, **kwargs):
        context = super(UserLibraryView, self).get_context_data(**kwargs)
        return context

    def get_queryset(self):
        papers = self.request.user.lib.papers.all()
        return papers

library = UserLibraryView.as_view()


class AjaxableResponseMixin(object):
    """
    Mixin to add AJAX support to a form.
    Must be used with an object-based FormView (e.g. UpdateView)
    """

    def form_invalid(self, form):
        response = super(AjaxableResponseMixin, self).form_invalid(form)
        if self.request.is_ajax():
            return JsonResponse(form.errors, status=400)
        else:
            return response

    def form_valid(self, form):
        # We make sure to call the parent's form_valid() method because
        # it might to do some processing (in the case of CreateView, it will
        # call form.save() for example)
        response = super(AjaxableResponseMixin, self).form_valid(form)
        if self.request.is_ajax():
            data = self.get_ajax_data()
            return JsonResponse(data)
        else:
            return response

    def get_ajax_data(self):
        raise NotImplementedError


class UserBasicInfoUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                              UpdateView):
    form_class = UserBasicNoEmailForm

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        return reverse('core:home')

    def get_ajax_data(self):
        data = {'first_name': self.request.user.first_name,
                'last_name': self.request.user.last_name}
        return data

update_basic_info = UserBasicInfoUpdateView.as_view()


class UserAffiliationUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                FormView):
    form_class = UserAffiliationForm

    def get_object(self, queryset=None):
        return self.request.user.affiliation

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        if Affiliation.objects.filter(**form.cleaned_data).exists():  # invalid because of uniqueness
            self.request.user.affiliation = \
                Affiliation.objects.get(**form.cleaned_data)
            self.request.user.save()
            return super(UserAffiliationUpdateView, self).form_valid(form)
        else:
            return super(UserAffiliationUpdateView, self).form_valid(form)

    def form_valid(self, form):
        affiliation = form.save()
        self.request.user.affiliation = affiliation
        self.request.user.save()
        return super(UserAffiliationUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'department': self.request.user.affiliation.department,
                'institution': self.request.user.affiliation.institution,
                'city': self.request.user.affiliation.city,
                'state': self.request.user.affiliation.state,
                'country': self.request.user.affiliation.country}
        return data

update_affiliation = UserAffiliationUpdateView.as_view()

