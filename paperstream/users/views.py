# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging


from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView, FormView, DetailView
from django.views.generic.edit import DeleteView
from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import JsonResponse
from django.utils import timezone
from django.template.loader import get_template
from django.core.mail import EmailMultiAlternatives
from django.template import Context

from braces.views import LoginRequiredMixin

from paperstream.core.views import BasePaperListView
from paperstream.core.mixins import AjaxableResponseMixin
from paperstream.library.models import Paper

from .forms import UserBasicForm, UserAffiliationForm, \
    UserAuthenticationForm, UserTrendSettingsForm, UserStreamSettingsForm, \
    UpdateUserNameForm, UpdateUserPositionForm, UpdateUserTitleForm, \
    UserEmailDigestSettingsForm
from .models import Affiliation, UserLibPaper, UserTaste, UserSettings
from .mixins import ProfileModalFormsMixin, SettingsModalFormsMixin
from .tasks import update_lib


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

    def get_ajax_data(self):
        data = {'email': self.request.user.email,
                'redirect': self.get_success_url(),
                }
        return data

ajax_signin = UserLoginView.as_view()


class UserBasicInfoSignupView(AjaxableResponseMixin, FormView):

    form_class = UserBasicForm
    template_name = 'user/signup_form.html'

    def get_success_url(self, **kwargs):
        return reverse('social:complete',
        kwargs={'backend': self.request.session['partial_pipeline']['backend']})

    def get_context_data(self, **kwargs):
        context = super(UserBasicInfoSignupView, self).get_context_data(**kwargs)
        context['stage'] = 'basic_info'
        context['has_email'] = True if self.initial.get('email') else False
        return context

    def get_initial(self):
        initial = super(UserBasicInfoSignupView, self).get_initial()
        details = self.request.session['partial_pipeline']['kwargs']['details']
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

    def get_ajax_data(self):
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


# User Library
# ---------------
class UserLibraryPaperListView(BasePaperListView):
    model = UserLibPaper
    # template_name = 'user/user_library.html'
    # page_template = 'user/user_library_sub_page.html'
    template_name = 'feeds/feed.html'
    page_template = 'feeds/feed_sub_page.html'

    def get_context_stats(self):
        context = super(UserLibraryPaperListView, self).get_context_stats()
        # Trash counter
        context['trash_counter'] = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib, is_trashed=True)\
            .count()
        # Like counter
        context['likes_counter'] = UserTaste.objects\
            .filter(user=self.request.user, is_liked=True)\
            .count()
        # library counter
        context['library_counter'] = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib, is_trashed=False)\
            .count()
        return context

    def get_context_data(self, **kwargs):
        context = super(UserLibraryPaperListView, self).get_context_data(**kwargs)
        context.update(self.get_context_endless(**kwargs))
        context.update(self.get_context_stats())
        context.update(self.get_context_usertaste())
        context.update(self.get_context_userlib())
        context.update(self.get_context_journals_filter())
        context.update(self.get_context_authors_filter())
        context.update(self.get_context_search_query())

        return context

    def get_queryset(self):

        # Retrieve get arguments if any
        self.parse_ajax_data()

        # Get data
        self.original_qs = self.get_original_queryset()

        # Exclude rejected paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')
        self.original_qs = self.original_qs.exclude(paper__in=papers_ticked)

        # filter
        queryset = self.original_qs
        if self.journals_filter:
            queryset = self.filter_journals(queryset)

        if self.authors_filter:
            queryset = self.filter_authors(queryset)

        if self.search_query:
            queryset = self.filter_search_query(queryset)

        return queryset

    def get_original_queryset(self):
        raise NotImplemented


class UserLibraryView(UserLibraryPaperListView):

    def get_original_queryset(self):
        return UserLibPaper.objects\
            .filter(userlib=self.request.user.lib,
                    is_trashed=False)

library = UserLibraryView.as_view()


class UserLibraryViewXML(UserLibraryView):
    # page_template = 'user/user_library_sub_page2.html'
    page_template = 'feeds/feed_sub_page2.html'

    def get_context_data(self, **kwargs):
        context = super(UserLibraryViewXML, self).get_context_data(**kwargs)
        context.update(self.get_context_filter_json(context))
        return context

library_xml = UserLibraryViewXML.as_view()


class UserLibraryTrashView(UserLibraryPaperListView):

    def get_original_queryset(self):
        return UserLibPaper.objects\
            .filter(userlib=self.request.user.lib,
                    is_trashed=True)

library_trash = UserLibraryTrashView.as_view()


class UserLibraryTrashViewXML(UserLibraryTrashView):
    # page_template = 'user/user_library_sub_page2.html'
    page_template = 'feeds/feed_sub_page2.html'

    def get_context_data(self, **kwargs):
        context = super(UserLibraryTrashViewXML, self).get_context_data(**kwargs)
        context.update(self.get_context_filter_json(context))
        return context

library_trash_xml = UserLibraryTrashViewXML.as_view()


class UserLibraryLikesView(UserLibraryPaperListView):

    def get_original_queryset(self):
        return UserTaste.objects\
            .filter(user=self.request.user,
                    is_liked=True)

library_likes = UserLibraryLikesView.as_view()


class UserLibraryLikesViewXML(UserLibraryLikesView):
    # page_template = 'user/user_library_sub_page2.html'
    page_template = 'feeds/feed_sub_page2.html'

    def get_context_data(self, **kwargs):
        context = super(UserLibraryLikesViewXML, self).get_context_data(**kwargs)
        context.update(self.get_context_filter_json(context))
        return context

library_likes_xml = UserLibraryLikesViewXML.as_view()


# Profile
# -------------------
class ProfileView(LoginRequiredMixin, ProfileModalFormsMixin, DetailView):
    model = User
    template_name = 'user/profile.html'

    def get_object(self, **kwargs):
        return self.request.user

    def get_context_data(self, **kwargs):
        context = super(ProfileView, self).get_context_data(**kwargs)
        context['library_counter'] = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib, is_trashed=False)\
            .count()
        context['likes_counter'] = UserTaste.objects\
            .filter(user=self.request.user, is_liked=True)\
            .count()
        return context

profile = ProfileView.as_view()


class SettingsView(LoginRequiredMixin, SettingsModalFormsMixin, DetailView):

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

    def get_ajax_data(self):
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

    def get_ajax_data(self):
        data = {'title': self.request.user.title}
        return data

update_title = UpdateUserTitleView.as_view()


class UpdateUserPositionView(LoginRequiredMixin, AjaxableResponseMixin, UpdateView):
    form_class = UpdateUserPositionForm

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        return reverse('core:home')

    def get_ajax_data(self):
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

    def get_ajax_data(self):
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


class UserStreamSettingsUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                   FormView):
    form_class = UserStreamSettingsForm

    def get_object(self, queryset=None):
        return self.request.user.settings

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        return super(UserStreamSettingsUpdateView, self).form_invalid(form)

    def form_valid(self, form):
        self.request.user.settings.stream_vector_weight = form.cleaned_data['stream_vector_weight']
        self.request.user.settings.stream_author_weight = form.cleaned_data['stream_author_weight']
        self.request.user.settings.stream_journal_weight = form.cleaned_data['stream_journal_weight']
        self.request.user.settings.save()
        return super(UserStreamSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'stream_vector_weight': '{0:.2f}'.format(self.request.user.settings.stream_vector_weight),
                'stream_author_weight': '{0:.2f}'.format(self.request.user.settings.stream_author_weight),
                'stream_journal_weight': '{0:.2f}'.format(self.request.user.settings.stream_journal_weight),
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

    def form_valid(self, form):
        self.request.user.settings.trend_doc_weight = form.cleaned_data['trend_doc_weight']
        self.request.user.settings.trend_altmetric_weight = form.cleaned_data['trend_altmetric_weight']
        self.request.user.settings.save()
        return super(UserTrendSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'trend_doc_weight': '{0:.2f}'.format(self.request.user.settings.trend_doc_weight),
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

    def get_ajax_data(self):
        data = {'email_digest_frequency': self.request.user.settings.get_email_digest_frequency_display(),
                }
        return data

update_email_digest_settings = UserEmailDigestSettingsUpdateView.as_view()


@login_required
def user_init_step(request):
    if request.method == 'GET':
        messages = []
        done = False
        redirect = ''
        if request.user.init_step == 'LIB':
            messages = ['Syncing your library with PubStream' ,
                        '({0})'.format(request.user.lib.count_papers)]
        elif request.user.init_step == 'STR':
            messages = ['Building your personalized streams',
                        '(1/2)']
        elif request.user.init_step == 'TRE':
            messages = ['Building your personalized streams',
                        '2/2']
        elif request.user.init_step == 'IDL':
            done = True
            messages = ['Done', '']
            redirect = reverse('feeds:stream')

        data = {'done': done,
                'step': request.user.init_step,
                'messages': messages,
                'redirect': redirect}

        return JsonResponse(data)

@login_required
def user_update_step(request):
    if request.method == 'GET':
        messages = []
        done = False
        redirect = ''
        if request.user.streams.first().state == 'ING':
            messages = ['Updating Your News']
        elif request.user.trends.first().state == 'ING':
            messages = ['Updating Top News']
        elif request.user.streams.first().state == 'IDL' and \
                        request.user.trends.first().state == 'IDL':
            done = True

        data = {'done': done,
                'messages': messages}

        return JsonResponse(data)


@login_required
def update_library(request):
    if request.is_ajax():
        if request.user.lib.state == 'IDL':
            provider_name = request.user.social_auth.first().provider
            update_lib.delay(request.user.pk, provider_name)
            data = {}
            return JsonResponse(data)
        else:
            data = {'errors': 'Error: Your library is already syncing'}
            return JsonResponse(data, status=409)


@login_required
def pin_call(request):
    if request.method == 'POST':
        pk = int(request.POST.get('pk'))
        source = request.POST.get('source', '')
        paper_ = get_object_or_404(Paper, pk=pk)
        if 'stream' in source:
            context_source = 'stream'
        elif 'trend' in source:
            context_source = 'trend'
        elif 'library' in source:
            context_source = 'library'
        else:
            context_source = source
        ut, _ = UserTaste.objects.get_or_create(
            paper=paper_,
            user=request.user)
        if ut.is_liked:
            ut.is_liked = False
        else:
            ut.is_liked = True
        ut.context_source = context_source
        ut.save()
        data = {'is_liked': ut.is_liked,
                'is_ticked': ut.is_ticked,
                'likes_counter': UserTaste.objects
                    .filter(user=request.user, is_liked=True)
                    .count()}
        return JsonResponse(data)
    else:
        return redirect('feeds:main')


@login_required
def ban_call(request):
    if request.method == 'POST':
        pk = int(request.POST.get('pk'))
        source = request.POST.get('source')
        paper_ = get_object_or_404(Paper, pk=pk)
        if 'stream' in source:
            context_source = 'stream'
        elif 'trend' in source:
            context_source = 'trend'
        elif 'library' in source:
            context_source = 'library'
        else:
            context_source = source
        ut, _ = UserTaste.objects.get_or_create(
            paper=paper_,
            user=request.user,
            context_source=context_source)
        if ut.is_ticked:
            ut.is_ticked = False
        else:
            ut.is_ticked = True
        ut.save()
        data = {'is_liked': ut.is_liked,
                'is_ticked': ut.is_ticked,
                'likes_counter': UserTaste.objects
                    .filter(user=request.user, is_liked=True)
                    .count()}
        return JsonResponse(data)
    else:
        return redirect('feeds:main')


@login_required
def add_call(request):
    """Push paper to user reference manager and add to local library"""
    if request.method == 'POST':
        pk = int(request.POST.get('pk'))
        user = request.user
        provider_name = user.social_auth.first().provider

        # get social
        social = user.social_auth.get(provider=provider_name)

        # get backend
        backend = social.get_backend_instance()

        # build session
        session = backend.get_session(social, user)

        # Get paper
        paper = Paper.objects.get(pk=pk)

        # push paper to lib
        err, paper_provider_id = backend.add_paper(session, paper)

        # return JSON data
        if not err:
            # add paper locally
            backend.associate_paper(paper, user, {'created': timezone.now().date()},
                                paper_provider_id)
            backend.associate_journal(paper.journal, user)
            data = {'success': True,
                    'trash_counter':  UserLibPaper.objects\
                        .filter(userlib=request.user.lib, is_trashed=True)\
                        .count(),
                    'library_counter': UserLibPaper.objects\
                        .filter(userlib=request.user.lib, is_trashed=False)\
                        .count(),
                    'likes_counter': UserTaste.objects\
                        .filter(user=request.user, is_liked=True)\
                        .count(),
                    'message': ''}
        else:
            data = {'success': False,
                    'message': 'Cannot add this paper to your library. Something went wrong'}
        return JsonResponse(data)


@login_required
def trash_call(request):
    """Trash paper from user library"""
    if request.method == 'POST':
        pk = int(request.POST.get('pk'))
        user = request.user

        # get social
        social = user.social_auth.first()

        # get backend
        backend = social.get_backend_instance()

        # build session
        session = backend.get_session(social, user)

        # Get userlibpaper
        ulp = user.lib.userlib_paper.get(paper_id=pk)

        # remove paper from provider lib
        err = backend.trash_paper(session, ulp)

        if not err:
            # remove paper locally from user library
            ulp.is_trashed = True
            ulp.save(update_fields=['is_trashed'])
            # remove from UserTaste
            ut, _ = UserTaste.objects.get_or_create(user=request.user,
                                                    paper_id=pk)
            ut.is_liked = False
            ut.is_ticked = True
            ut.save()
            # build json data
            data = {'success': True,
                    'trash_counter':  UserLibPaper.objects\
                        .filter(userlib=request.user.lib, is_trashed=True)\
                        .count(),
                    'library_counter': UserLibPaper.objects\
                        .filter(userlib=request.user.lib, is_trashed=False)\
                        .count(),
                    'likes_counter': UserTaste.objects\
                        .filter(user=request.user, is_liked=True)\
                        .count(),
                    'message': ''}
        else:
            data = {'success': False,
                    'message': 'Cannot remove this paper to your library. Something went wrong'}
            logger.error('Fail trashing paper {pk} for user {pk_user}'.format(
                pk=pk,
                pk_user=user.pk))
        return JsonResponse(data)


@login_required
def restore_call(request):
    """Restore paper from trash"""
    if request.method == 'POST':
        pk = int(request.POST.get('pk'))
        user = request.user
        provider_name = user.social_auth.first().provider

        # get social
        social = user.social_auth.get(provider=provider_name)

        # get backend
        backend = social.get_backend_instance()

        # build session
        session = backend.get_session(social, user)

        # Get paper
        paper = Paper.objects.get(pk=pk)

        # push paper to lib
        err, paper_provider_id = backend.add_paper(session, paper)

        # return JSON data
        if not err:
            # restore paper locally from user library
            backend.associate_paper(paper, user, {'created': timezone.now().date()},
                                    paper_provider_id)
            backend.associate_journal(paper.journal, user)
            data = {'success': True,
                    'trash_counter':  UserLibPaper.objects\
                        .filter(userlib=request.user.lib, is_trashed=True)\
                        .count(),
                    'library_counter': UserLibPaper.objects\
                        .filter(userlib=request.user.lib, is_trashed=False)\
                        .count(),
                    'likes_counter': UserTaste.objects\
                        .filter(user=request.user, is_liked=True)\
                        .count(),
                    'message': ''}
        else:
            data = {'success': False,
                    'message': 'Cannot add this paper to your library. Something went wrong'}
        return JsonResponse(data)


@login_required
def send_invite(request):
    if request.POST:
        email_to = request.POST.get('email')

        subject = 'An invitation to try Etalia'
        to = [email_to]
        from_email = request.user.email
        ctx = {'bucket_url': settings.EMAIL_STATIC_BUCKET,
               'root_url': request.META.get('HTTP_ORIGIN')}
        text_content = ''
        html_content = get_template(settings.INVITE_EMAIL_TEMPLATE)\
            .render(Context(ctx))
        email = EmailMultiAlternatives(subject, text_content, to=to, from_email=from_email)
        email.attach_alternative(html_content, "text/html")
        # email.content_subtype = 'html'
        email.send()

        return JsonResponse(data={'success': True})
    else:
        redirect('invite:home')