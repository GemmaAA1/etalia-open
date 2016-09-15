# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging

from celery.canvas import chain

from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView, FormView, DetailView, RedirectView
from django.views.generic.edit import DeleteView
from django.conf import settings
from django.contrib.auth import get_user_model
from django.utils.text import slugify
from django.http import JsonResponse
from django.db import transaction
from django.views.decorators.cache import never_cache

from braces.views import LoginRequiredMixin

from etalia.core.mixins import AjaxableResponseMixin, NavFlapMixin
from etalia.nlp.tasks import update_userfingerprint
from etalia.feeds.tasks import update_trend, update_stream, update_threadfeed
from .forms import UserBasicForm, UserAffiliationForm, \
    UserAuthenticationForm, UserTrendSettingsForm, UserStreamSettingsForm, \
    UpdateUserNameForm, UpdateUserPositionForm, UpdateUserTitleForm, \
    UserEmailDigestSettingsForm, UserFingerprintSettingsForm
from etalia.library.models import PaperUser
from etalia.library.constants import PAPER_PINNED
from .models import Affiliation, UserLibPaper, UserSettings
from .mixins import ProfileModalFormsMixin, SettingsModalFormsMixin
from .tasks import update_lib
from .constants import INIT_STEPS, USERLIB_IDLE, USERLIB_SYNCING
from .utils import send_invite_email

logger = logging.getLogger(__name__)
User = get_user_model()


def logout(request):
    """Logs out user"""
    auth_logout(request)
    return redirect('/')


def done(request):
    redirect('landing')


# User authentication
# ---------------
class UserLoginView(AjaxableResponseMixin, FormView):

    form_class = UserAuthenticationForm
    redirect_field_name = settings.LOGIN_REDIRECT_URL

    def form_valid(self, form):
        login(self.request, form.get_user())
        response = super(AjaxableResponseMixin, self).form_valid(form)
        if self.request.is_ajax():
            data = self.get_ajax_data()
            return JsonResponse(data)
        else:
            return response

    def get_success_url(self):
        return reverse('feeds:main')

    def get_ajax_data(self, *args, **kwargs):
        data = {'email': self.request.user.email,
                'redirect': self.get_success_url(),
                }
        return data

ajax_signin = UserLoginView.as_view()


class UserBasicInfoSignupView(AjaxableResponseMixin, FormView):

    form_class = UserBasicForm
    template_name = 'user/signup_form.html'

    def get_success_url(self, **kwargs):
        return reverse(
            'social:complete',
            kwargs={
                'backend': self.request.session['partial_pipeline']['backend']
            })

    def get_context_data(self, **kwargs):
        context = super(UserBasicInfoSignupView, self).get_context_data(**kwargs)
        context['stage'] = 'basic_info'
        context['has_email'] = True if self.initial.get('email') else False
        return context

    def get_initial(self):
        initial = super(UserBasicInfoSignupView, self).get_initial()
        try:
            details = self.request.session['partial_pipeline']['kwargs']['details']
        except KeyError:
            details = {'first_name': '', 'last_name': '', 'email': ''}
        initial['first_name'] = details.get('first_name', '')
        initial['last_name'] = details.get('last_name', '')
        initial['email'] = details.get('email', '')
        self.initial = initial
        return initial

    def form_valid(self, form):
        self.request.session['basic_info'] = {
            'first_name': form.cleaned_data['first_name'],
            'last_name': form.cleaned_data['last_name'],
            'email': form.cleaned_data['email'],
            }
        return super(UserBasicInfoSignupView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'redirect': self.get_success_url()}
        return data

require_basic_info = UserBasicInfoSignupView.as_view()


class UserAffiliationSignupView(FormView):

    form_class = UserAffiliationForm
    template_name = 'user/signup_form.html'

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
        return dict(details.get('tmp_affiliation', ''), **initial)

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
    return render(request, 'user/signup_form.html', context)


# -----------------------------------------------------------------------------
#  USER PROFILE
# -----------------------------------------------------------------------------


class ProfileView(LoginRequiredMixin, ProfileModalFormsMixin, NavFlapMixin,
                  DetailView):
    model = User
    template_name = 'user/profile.html'

    def get_object(self, **kwargs):
        return self.request.user

    def get_context_data(self, **kwargs):
        context = super(ProfileView, self).get_context_data(**kwargs)
        context['library_counter'] = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib)\
            .count()
        context['likes_counter'] = PaperUser.objects\
            .filter(user=self.request.user, watch=PAPER_PINNED)\
            .count()
        return context

profile = ProfileView.as_view()


class UserProfileSlugView(NavFlapMixin, DetailView):

    model = User
    template_name = 'user/user_profile.html'

user_profile_slug = UserProfileSlugView.as_view()


class UserProfilePkView(RedirectView):
    """Redirect to slug paper url"""

    permanent = False
    query_string = True
    pattern_name = 'user:user-profile-slug'

    def get_redirect_url(self, *args, **kwargs):
        user = User.objects.get(pk=kwargs['pk'])
        kwargs['slug'] = slugify(' '.join([user.first_name, user.last_name]))
        return super(UserProfilePkView, self).get_redirect_url(*args, **kwargs)

user_profile_pk = UserProfilePkView.as_view()


class SettingsView(LoginRequiredMixin, SettingsModalFormsMixin, NavFlapMixin,
                   DetailView):

    model = UserSettings
    template_name = 'user/settings.html'

    def get_object(self, **kwargs):
        return self.request.user.settings

settings_view = SettingsView.as_view()


class UpdateUserNameView(LoginRequiredMixin, AjaxableResponseMixin, UpdateView):
    form_class = UpdateUserNameForm

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        return reverse('core:home')

    def get_ajax_data(self, *args, **kwargs):
        data = {'first_name': self.request.user.first_name,
                'last_name': self.request.user.last_name}
        return data

update_name = UpdateUserNameView.as_view()


class UpdateUserTitleView(LoginRequiredMixin, AjaxableResponseMixin, UpdateView):
    form_class = UpdateUserTitleForm

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        return reverse('core:home')

    def get_ajax_data(self, *args, **kwargs):
        data = {'title': self.request.user.title}
        return data

update_title = UpdateUserTitleView.as_view()


class UpdateUserPositionView(LoginRequiredMixin, AjaxableResponseMixin, UpdateView):
    form_class = UpdateUserPositionForm

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        return reverse('core:home')

    def get_ajax_data(self, *args, **kwargs):
        data = {'position': self.request.user.position}
        return data

update_position = UpdateUserPositionView.as_view()


class UserAffiliationUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                FormView):
    form_class = UserAffiliationForm

    def get_object(self, queryset=None):
        return self.request.user.affiliation

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        if Affiliation.objects.filter(**form.cleaned_data).exists():
            self.request.user.affiliation = \
                Affiliation.objects.get(**form.cleaned_data)
            self.request.user.save()
            return super(UserAffiliationUpdateView, self).form_valid(form)
        else:
            return super(UserAffiliationUpdateView, self).form_invalid(form)

    def form_valid(self, form):
        affiliation = form.save()
        self.request.user.affiliation = affiliation
        self.request.user.save()
        return super(UserAffiliationUpdateView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'print-affiliation':
                    self.request.user.affiliation.print_affiliation}
        return data

update_affiliation = UserAffiliationUpdateView.as_view()


class UserDeleteView(LoginRequiredMixin, ProfileModalFormsMixin, DeleteView):
    model = User
    template_name = 'user/user_confirm_delete.html'
    success_url = reverse_lazy('core:home')

    def get_object(self, **kwargs):
        return self.request.user

delete_user = UserDeleteView.as_view()


class UserFingerprintSettingsUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                        FormView):
    form_class = UserFingerprintSettingsForm

    def get_object(self, queryset=None):
        return self.request.user.settings

    def get_success_url(self):
        return reverse('core:home')

    @transaction.non_atomic_requests
    def form_valid(self, form):
        self.request.user.settings.fingerprint_roll_back_deltatime = \
            form.cleaned_data['fingerprint_roll_back_deltatime']
        self.request.user.settings.save()
        # Trigger update
        task = chain(
            update_userfingerprint.s(self.request.user.id),
            update_stream.s(),
            update_trend.s(),
            update_threadfeed.s(),
        )
        task()
        return super(UserFingerprintSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'fingerprint_roll_back_deltatime': '{0:d}'.format(self.request.user.settings.fingerprint_roll_back_deltatime),
                }
        return data

update_fingerprint_settings = UserFingerprintSettingsUpdateView.as_view()


class UserStreamSettingsUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                   FormView):
    form_class = UserStreamSettingsForm

    def get_object(self, queryset=None):
        return self.request.user.settings

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        return super(UserStreamSettingsUpdateView, self).form_invalid(form)

    @transaction.non_atomic_requests
    def form_valid(self, form):
        self.request.user.settings.stream_vector_weight = form.cleaned_data['stream_vector_weight']
        self.request.user.settings.stream_author_weight = form.cleaned_data['stream_author_weight']
        self.request.user.settings.stream_journal_weight = form.cleaned_data['stream_journal_weight']
        self.request.user.settings.stream_score_threshold = form.cleaned_data['stream_score_threshold']
        self.request.user.settings.save()
        update_stream.delay(self.request.user.id)
        return super(UserStreamSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'stream_vector_weight': '{0:.2f}'.format(self.request.user.settings.stream_vector_weight),
                'stream_author_weight': '{0:.2f}'.format(self.request.user.settings.stream_author_weight),
                'stream_journal_weight': '{0:.2f}'.format(self.request.user.settings.stream_journal_weight),
                'stream_score_threshold': '{0:.2f}'.format(self.request.user.settings.stream_score_threshold),
                }
        return data

update_stream_settings = UserStreamSettingsUpdateView.as_view()


class UserTrendSettingsUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                  FormView):
    form_class = UserTrendSettingsForm

    def get_object(self, queryset=None):
        return self.request.user.settings

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        return super(UserTrendSettingsUpdateView, self).form_invalid(form)

    @transaction.non_atomic_requests
    def form_valid(self, form):
        self.request.user.settings.trend_doc_weight = form.cleaned_data['trend_doc_weight']
        self.request.user.settings.trend_altmetric_weight = form.cleaned_data['trend_altmetric_weight']
        self.request.user.settings.trend_score_threshold = form.cleaned_data['trend_score_threshold']
        self.request.user.settings.save()
        update_trend.delay(self.request.user.id)
        return super(UserTrendSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'trend_score_threshold': '{0:.2f}'.format(self.request.user.settings.trend_score_threshold),
                'trend_doc_weight': '{0:.2f}'.format(self.request.user.settings.trend_doc_weight),
                'trend_altmetric_weight': '{0:.2f}'.format(self.request.user.settings.trend_altmetric_weight),
                }
        return data

update_trend_settings = UserTrendSettingsUpdateView.as_view()


class UserEmailDigestSettingsUpdateView(LoginRequiredMixin,
                                        AjaxableResponseMixin,
                                        FormView):
    form_class = UserEmailDigestSettingsForm

    def get_object(self, queryset=None):
        return self.request.user.settings

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        return super(UserEmailDigestSettingsUpdateView, self).form_invalid(form)

    def form_valid(self, form):
        self.request.user.settings.email_digest_frequency = form.cleaned_data['email_digest_frequency']
        self.request.user.settings.save()
        return super(UserEmailDigestSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'email_digest_frequency': self.request.user.settings.get_email_digest_frequency_display(),
                }
        return data

update_email_digest_settings = UserEmailDigestSettingsUpdateView.as_view()


@login_required
@never_cache
def user_init_check(request):
    if request.method == 'GET':
        messages = []
        done = False
        steps = dict(INIT_STEPS)
        redirect = ''
        for k, v in steps.items():
            if request.user.init_step == k:
                messages = [v, '']

                # special case
                if request.user.init_step == 'LIB':
                    messages[1] = '({0} papers)'.format(request.user.lib.count_papers)

                if request.user.init_step == 'IDL':
                    done = True
                    messages = ['Done', '']
                    redirect = reverse('feeds:my_feeds')

        data = {'done': done,
                'step': request.user.init_step,
                'messages': messages,
                'redirect': redirect}

        return JsonResponse(data)


@login_required
def user_update_fingerprint_check(request):
    if request.method == 'GET':
        if request.user.fingerprint.first().state == 'ING':
            messages = ['Updating Your Fingerprint', '']
            done = False
        else:
            messages = []
            done = True
        data = {'done': done,
                'messages': messages}
        return JsonResponse(data)


@login_required
def user_update_stream_check(request):
    if request.method == 'GET':
        if request.user.streams.first().state == 'ING':
            messages = ['Updating Your Stream', '']
            done = False
        else:
            messages = []
            done = True
        data = {'done': done,
                'messages': messages}
        return JsonResponse(data)


@login_required
def user_update_trend_check(request):
    if request.method == 'GET':
        if request.user.trends.first().state == 'ING':
            messages = ['Updating Trends', '']
            done = False
        else:
            messages = []
            done = True

        data = {'done': done,
                'messages': messages}
        return JsonResponse(data)


@login_required
def user_update_settings_check(request):
    if request.method == 'GET':
        if request.user.streams.first().state == 'ING':
            return user_update_stream_check(request)
        elif request.user.trend.first().state == 'ING':
            return user_update_trend_check(request)
        elif request.user.lib.state == USERLIB_SYNCING:
            return user_update_library_check(request)
        elif request.user.fingerprint.first().state == 'ING':
            return user_update_fingerprint_check(request)
        else:
            data = {'done': True, 'messages': []}
            return JsonResponse(data)


@login_required
def user_update_library_check(request):
    if request.method == 'GET':
        if request.user.lib.state == USERLIB_SYNCING:
            messages = ['Updating You Library', '']
            done = False
        else:
            messages = []
            done = True

        data = {'done': done,
                'messages': messages}
        return JsonResponse(data)


@login_required
def update_library(request):
    if request.is_ajax():
        if request.user.lib.state == USERLIB_IDLE:
            provider_name = request.user.social_auth.first().provider
            update_lib.delay(request.user.pk, provider_name)
            data = {}
            return JsonResponse(data)
        else:
            data = {'errors': 'Error: Your library is already syncing'}
            return JsonResponse(data, status=409)


@login_required
def send_invite(request):
    if request.POST:
        email_to = request.POST.get('email')
        send_invite_email(
            email_to=email_to,
            on_behalf=request.user,
        )
        return JsonResponse(data={'success': True})
    else:
        redirect('invite:home')


def test_template(request):
    from django.template.response import TemplateResponse
    papers = request.user.streams.first().papers.all()
    papers = papers.select_related('altmetric')
    papers = papers.select_related('journal')
    papers = papers.prefetch_related('authors')
    papers = list(papers[:settings.PERIODIC_RECOMMENDATION_NUMBER_PAPERS])

    return TemplateResponse(
        request,
        settings.PERIODIC_RECOMMENDATION_TEMPLATE,
        {'papers': papers,
         'bucket_url': settings.EMAIL_STATIC_BUCKET,
         'root_url': 'http://etalia.io',
         }
    )





