# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import logging


from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView, FormView, DetailView
from django.views.generic.edit import DeleteView
from django.views.generic.base import TemplateView
from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import JsonResponse
from django.utils import timezone
from django.template.loader import get_template
from django.core.mail import EmailMultiAlternatives
from django.template import Context

from braces.views import LoginRequiredMixin

from etalia.core.views import BasePaperListView
from etalia.core.mixins import AjaxableResponseMixin, NavFlapMixin, \
    XMLMixin

from .forms import UserBasicForm, UserAffiliationForm, \
    UserAuthenticationForm, UserTrendSettingsForm, UserStreamSettingsForm, \
    UpdateUserNameForm, UpdateUserPositionForm, UpdateUserTitleForm, \
    UserEmailDigestSettingsForm, UserTasteForm, UserLibPaperForm
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
#  USER LIBRARY
# -----------------------------------------------------------------------------


class UserLibraryPaperListView(BasePaperListView):
    model = UserLibPaper
    template_name = 'user/user_library.html'
    page_template = 'user/user_library_sub_page.html'

    def get_context_stats(self):
        context = super(UserLibraryPaperListView, self).get_context_stats()
        # Trash counter
        context['trash_counter'] = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib, is_trashed=True)\
            .count()
        # Like counter
        context['likes_counter'] = UserTaste.objects\
            .filter(user=self.request.user, is_pinned=True)\
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
        context.update(self.get_context_control_session())

        return context

    def get_queryset(self):

        # Retrieve get arguments if any
        self.get_input_data()

        # Get data
        self.original_qs = self.get_original_queryset()

        # Exclude rejected paper
        papers_banned = UserTaste.objects\
            .filter(user=self.request.user,
                    is_banned=True)\
            .values('paper')
        self.original_qs = self.original_qs.exclude(paper__in=papers_banned)

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


class BaseUserLibraryView(UserLibraryPaperListView):

    control_session = 'control_library'

    def get_original_queryset(self):
        return UserLibPaper.objects\
            .filter(userlib=self.request.user.lib,
                    is_trashed=False)\
            .select_related('paper',
                            'paper__journal',
                            'paper__altmetric')

    def get_context_data(self, **kwargs):
        context = super(BaseUserLibraryView, self).get_context_data(**kwargs)
        context['tab'] = 'main'
        return context


class UserLibraryView(LoginRequiredMixin, NavFlapMixin, BaseUserLibraryView):
    pass

library = UserLibraryView.as_view()


class UserLibraryViewXML(LoginRequiredMixin, NavFlapMixin, XMLMixin,
                         BaseUserLibraryView):
    page_template = 'user/user_library_sub_page2.html'

library_xml = UserLibraryViewXML.as_view()


class BaseUserLibraryTrashView(UserLibraryPaperListView):

    control_session = 'control_trash'

    def get_original_queryset(self):
        return UserLibPaper.objects\
            .filter(userlib=self.request.user.lib,
                    is_trashed=True)\
            .select_related('paper',
                            'paper__journal',
                            'paper__altmetric')

    def get_context_data(self, **kwargs):
        context = super(BaseUserLibraryTrashView, self).get_context_data(**kwargs)
        context['tab'] = 'trash'
        return context


class UserLibraryTrashView(LoginRequiredMixin, NavFlapMixin,
                           BaseUserLibraryTrashView):
    pass

library_trash = UserLibraryTrashView.as_view()


class UserLibraryTrashViewXML(LoginRequiredMixin, NavFlapMixin, XMLMixin,
                              BaseUserLibraryTrashView):
    page_template = 'user/user_library_sub_page2.html'

library_trash_xml = UserLibraryTrashViewXML.as_view()


class BaseUserLibraryPinsView(UserLibraryPaperListView):

    control_session = 'control_pins'

    def get_original_queryset(self):
        return UserTaste.objects\
            .filter(user=self.request.user,
                    is_pinned=True)\
            .select_related('paper',
                            'paper__journal',
                            'paper__altmetric')

    def get_context_data(self, **kwargs):
        context = super(BaseUserLibraryPinsView, self).get_context_data(**kwargs)
        context['tab'] = 'pin'
        return context


class UserLibraryPinsView(LoginRequiredMixin, NavFlapMixin,
                          BaseUserLibraryPinsView):
    pass

library_pins = UserLibraryPinsView.as_view()


class UserLibraryPinsViewXML(LoginRequiredMixin, NavFlapMixin, XMLMixin,
                             BaseUserLibraryPinsView):
    page_template = 'user/user_library_sub_page2.html'

library_pins_xml = UserLibraryPinsViewXML.as_view()

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
            .filter(userlib=self.request.user.lib, is_trashed=False)\
            .count()
        context['likes_counter'] = UserTaste.objects\
            .filter(user=self.request.user, is_pinned=True)\
            .count()
        return context

profile = ProfileView.as_view()


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
        self.request.user.settings.stream_roll_back_deltatime = form.cleaned_data['stream_roll_back_deltatime']
        self.request.user.settings.save()
        return super(UserStreamSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        data = {'stream_vector_weight': '{0:.2f}'.format(self.request.user.settings.stream_vector_weight),
                'stream_author_weight': '{0:.2f}'.format(self.request.user.settings.stream_author_weight),
                'stream_journal_weight': '{0:.2f}'.format(self.request.user.settings.stream_journal_weight),
                'stream_roll_back_deltatime': '{0:d} months'.format(self.request.user.settings.stream_roll_back_deltatime),
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

    def get_ajax_data(self, *args, **kwargs):
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

    def get_ajax_data(self, *args, **kwargs):
        data = {'email_digest_frequency': self.request.user.settings.get_email_digest_frequency_display(),
                }
        return data

update_email_digest_settings = UserEmailDigestSettingsUpdateView.as_view()


@login_required
def user_init_check(request):
    if request.method == 'GET':
        messages = []
        done = False
        redirect = ''
        if request.user.init_step == 'LIB':
            messages = ['Syncing your library',
                        '({0} papers)'.format(request.user.lib.count_papers)]
        elif request.user.init_step == 'STR':
            messages = ['Building your streams',
                        '(1/2)']
        elif request.user.init_step == 'TRE':
            messages = ['Building your streams',
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
        elif request.user.streams.first().state == 'ING':
            return user_update_stream_check(request)
        elif request.user.lib.state == 'ING':
            return user_update_library_check(request)
        else:
            data = {'done': True, 'messages': []}
            return JsonResponse(data)

@login_required
def user_update_library_check(request):
    if request.method == 'GET':
        if request.user.lib.state == 'ING':
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
        if request.user.lib.state == 'IDL':
            provider_name = request.user.social_auth.first().provider
            update_lib.delay(request.user.pk, provider_name)
            data = {}
            return JsonResponse(data)
        else:
            data = {'errors': 'Error: Your library is already syncing'}
            return JsonResponse(data, status=409)


class UserPaperCallView(LoginRequiredMixin, AjaxableResponseMixin, FormView):
    """Abstract class to deal with intercation between user and paper"""
    def get_success_url(self):
        return reverse('core:home')

    def get_form_kwargs(self):
        kwargs = super(UserPaperCallView, self).get_form_kwargs()
        return kwargs

    def get_ajax_data(self, *args, **kwargs):
        return {
           'state': self.request.user.get_paper_state(
               kwargs['form'].cleaned_data['paper'].id),
           'counter': self.request.user.get_counters()
        }

    def get_session_backend(self):
        user = self.request.user
        provider_name = user.social_auth.first().provider
        # get social
        social = user.social_auth.get(provider=provider_name)
        # get backend
        backend = social.get_backend_instance()
        # build session
        session = backend.get_session(social, user)
        return session, backend


class PinCallView(UserPaperCallView):

    form_class = UserTasteForm

    def form_valid(self, form):
        ut, _ = UserTaste.objects.get_or_create(
            paper=form.cleaned_data['paper'],
            user=self.request.user
        )
        ut.is_pinned = not ut.is_pinned
        ut.source = form.cleaned_data['source']
        ut.save()
        return super(PinCallView, self).form_valid(form)

pin_call = PinCallView.as_view()


class BanCallView(UserPaperCallView):

    form_class = UserTasteForm

    def form_valid(self, form):
        ut, _ = UserTaste.objects.get_or_create(
            paper=form.cleaned_data['paper'],
            user=self.request.user
        )
        ut.is_banned = not ut.is_banned
        ut.source = form.cleaned_data['source']
        ut.save()
        return super(BanCallView, self).form_valid(form)

ban_call = BanCallView.as_view()


class AddCallView(UserPaperCallView):

    form_class = UserLibPaperForm

    def form_valid(self, form):
        session, backend = self.get_session_backend()
        # add paper to user lib
        paper = form.cleaned_data['paper']
        err, paper_provider_id = backend.add_paper(session, paper)
        if not err:
            backend.associate_paper(paper, self.request.user,
                                    {'created': timezone.now().date()},
                                    paper_provider_id)
            backend.associate_journal(paper.journal, self.request.user)
            return super(AddCallView, self).form_valid(form)
        else:
            data = {'success': False,
                    'message': 'Failed to add {title} to y{provider} '
                               'library'.format(title=paper.title,
                                                provider=backend.name)}
            return JsonResponse(data)

add_call = AddCallView.as_view()


class TrashCallView(UserPaperCallView):

    form_class = UserLibPaperForm

    def form_valid(self, form):
        session, backend = self.get_session_backend()
        # remove paper from provider lib
        ulp = self.request.user.lib.userlib_paper.get(
            paper=form.cleaned_data['paper'])
        err = backend.trash_paper(session, ulp)
        if not err:
            # remove paper locally from user library
            ulp.is_trashed = True
            ulp.save(update_fields=['is_trashed'])
            return super(TrashCallView, self).form_valid(form)
        else:
            data = {'success': False,
                    'message': 'Failed to trash {title} from your {provider} '
                               'library'.format(title=ulp.paper.title,
                                                provider=backend.name)}
            return JsonResponse(data)

trash_call = TrashCallView.as_view()


@login_required
def empty_trash_call(request):
    if request.method == 'POST':
        ulps = request.user.lib.userlib_paper.filter(is_trashed=True)
        ulps.delete()
        if request.is_ajax():
            data = {'counter': request.user.get_counters()}
            return JsonResponse(data)
        else:
            redirect('user:library')


class RestoreCallView(UserPaperCallView):

    form_class = UserLibPaperForm

    def form_valid(self, form):
        session, backend = self.get_session_backend()
        # trash paper to user lib
        paper = form.cleaned_data['paper']
        err, paper_provider_id = backend.add_paper(session, paper)
        if not err:
            # remove paper locally from user library
            ulp = self.request.user.lib.userlib_paper.get(paper=paper)
            ulp.is_trashed = False
            ulp.paper_provider_id = paper_provider_id
            ulp.save(update_fields=['is_trashed', 'paper_provider_id'])
            return super(RestoreCallView, self).form_valid(form)
        else:
            data = {'success': False,
                    'message': 'Failed to restore {title} to your {provider} '
                               'library'.format(title=paper.title,
                                                provider=backend.name)}
            return JsonResponse(data)

restore_call = RestoreCallView.as_view()

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


class ThreadsView(LoginRequiredMixin, NavFlapMixin, TemplateView):
    template_name = 'threads/threads.html'

threads = ThreadsView.as_view()
