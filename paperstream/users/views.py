# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import operator
import json
from functools import reduce
import logging
from collections import Counter


from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect, get_object_or_404
from django.db.models import Q
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView, FormView, DetailView
from django.views.generic.edit import DeleteView
from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import JsonResponse
from django.utils import timezone
from django.db.models.functions import Coalesce
from django.template.loader import render_to_string, get_template
from django.core.mail import EmailMessage
from django.core.mail import EmailMultiAlternatives
from django.template import Context

from braces.views import LoginRequiredMixin

from endless_pagination.views import AjaxListView

from paperstream.core.mixins import AjaxableResponseMixin, ModalMixin
from paperstream.library.models import Paper, Author

from .forms import UserBasicForm, UserAffiliationForm, \
    UserAuthenticationForm, UserTrendSettingsForm, UserStreamSettingsForm, \
    UpdateUserNameForm, UpdateUserPositionForm, UpdateUserTitleForm, \
    UserEmailDigestSettingsForm
from .models import Affiliation, UserLibPaper, UserTaste, FeedLayout, \
    UserSettings
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
class UserLibraryView(LoginRequiredMixin, ModalMixin, AjaxListView):
    model = UserLibPaper
    template_name = 'user/user_library.html'
    page_template = 'user/user_library_sub_page.html'
    first_page = 30
    per_page = 20
    context_object_name = 'ulp_list'
    list_n_ujl = 10
    size_max_journal_filter = 40
    size_load_journal_filter = 10
    size_max_author_filter = 40
    size_load_author_filter = 10
    journals_filter = None
    journals_filter_flag = 'all'
    authors_filter = None
    authors_filter_flag = 'all'
    sorting_flag = 'relevant'
    like_flag = False

    def update_from_filter(self):
        ufl, new = FeedLayout.objects.get_or_create(user=self.request.user)
        if new:
            ufl.stream_filter = {'journals_flag': self.journals_filter_flag,
                                 'authors_flag': self.authors_filter_flag,
                                 'sorting_flag': self.sorting_flag}
            ufl.trend_filter = {'journals_flag': self.journals_filter_flag,
                                 'authors_flag': self.authors_filter_flag,
                                 'sorting_flag': self.sorting_flag}
            ufl.library_filter = {'journals_flag': self.journals_filter_flag,
                                  'authors_flag': self.authors_filter_flag,
                                  'sorting_flag': self.sorting_flag}
            ufl.save()
        else:
            self.journals_filter = ufl.library_filter.get('journals')
            self.authors_filter = ufl.library_filter.get('authors')
            self.authors_filter_flag = ufl.library_filter.get('authors_flag') or self.authors_filter_flag
            self.journals_filter_flag = ufl.library_filter.get('journals_flag') or self.journals_filter_flag
            self.sorting_flag = ufl.library_filter.get('sorting_flag') or self.sorting_flag

        # From ajaxable filter
        if self.request.is_ajax():
            try:
                data = json.loads(list(self.request.GET)[0])
                if data['action'] == 'filter':
                    # get data
                    self.journals_filter = data.get('journals')
                    self.authors_filter = data.get('authors')
                    self.authors_filter_flag = data.get('authors_flag') or self.authors_filter_flag
                    self.journals_filter_flag = data.get('journals_flag') or self.journals_filter_flag
                    self.sorting_flag = data.get('sorting_flag') or self.sorting_flag
                    self.like_flag = data.get('like_flag')
                    ufl.library_filter = {
                        'journals': self.journals_filter,
                        'authors': self.authors_filter,
                        'journals_flag': self.journals_filter_flag,
                        'authors_flag': self.authors_filter_flag,
                        'sorting_flag': self.sorting_flag,
                    }
                    ufl.save()
            except ValueError:  # likely data from AjaxListView
                pass

    def get_context_data(self, **kwargs):
        context = super(AjaxListView, self).get_context_data(**kwargs)
        context = dict(context,
                       **super(UserLibraryView, self).get_context_data(**kwargs))
        context['first_page'] = self.first_page
        context['per_page'] = self.per_page
        if self.request.GET.get("query"):
            context['current_query'] = self.request.GET.get("query")

        # Get data for filter
        data = self.request.user.lib.papers.values('pk', 'journal__title', 'authors')
        j_titles = []
        authors = []
        check_papers = []
        for d in data:
            if d['journal__title'] and d['pk'] not in check_papers:  # rows are for different authors
                j_titles.append(d['journal__title'])
                check_papers.append(d['pk'])
            authors.append(d['authors'])

        # journal filter
        journals_counter = Counter(j_titles)
        j_titles = sorted(journals_counter, key=journals_counter.get,
                          reverse=True)[:self.size_max_journal_filter]
        # group authors by block identified with block tag
        block = 0
        context['journals'] = []
        for i, title in enumerate(j_titles):
            if not i % self.size_load_journal_filter:
                block += 1
            context['journals'].append((title, journals_counter[title], block))
        if len(context['journals']) <= self.size_load_journal_filter:
            context['journals_show_more'] = False
        else:
            context['journals_show_more'] = True

        # author filter
        authors_counter = Counter(authors)
        authors_pks = sorted(authors_counter, key=authors_counter.get,
                             reverse=True)
        clauses = ' '.join(['WHEN id=%s THEN %s' %
                            (pk, i) for i, pk in enumerate(authors_pks)])
        ordering = 'CASE %s END' % clauses

        authors_ordered = Author.objects\
            .filter(pk__in=authors_pks)\
            .exclude(Q(first_name='') & Q(last_name=''))\
            .extra(select={'ordering': ordering},
                   order_by=('ordering',))[:self.size_max_author_filter]
        # group authors by block identified with block tag
        block = 0
        context['authors'] = []
        for i, auth in enumerate(authors_ordered):
            if not i % 10:
                block += 1
            context['authors'].append((auth, auth.print_full, block))
        if len(context['authors']) <= self.size_load_author_filter:
            context['authors_show_more'] = False
        else:
            context['authors_show_more'] = True

        # Get user tastes label
        user_taste = UserTaste.objects\
            .filter(user=self.request.user)\
            .values_list('paper_id', 'is_liked')
        # reformat to dict
        user_taste = dict((key, {'liked': v1})
                          for key, v1 in user_taste)
        context['user_taste'] = user_taste

        # Set filter context
        if self.request.GET.get('query'):
            context['get_query'] = self.request.GET.get('query')

        if self.journals_filter:
            context['journals_filter'] = [title.lower() for title, val in self.journals_filter if val]
        else:
            context['journals_filter'] = [title.lower() for title in j_titles]
        if self.authors_filter:
            context['authors_filter'] = [pk for pk, val in self.authors_filter if val]
        else:
            context['authors_filter'] = authors_pks

        context['sorting_flag'] = self.sorting_flag

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

        # library tab
        context['tab'] = 'library'

        return context

    def get_queryset(self):
        query_set = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib,
                    is_trashed=False)
        return self.filter_queryset(query_set)

    def filter_queryset(self, queryset):
        # load stored filter
        self.update_from_filter()

        # journals
        if self.journals_filter_flag == 'all':
            if self.journals_filter:
                subset = []
                for journal, val in self.journals_filter:
                    if not val:
                        subset.append(Q(paper__journal__title__icontains=journal))
                if subset:
                    queryset = queryset.exclude(reduce(operator.or_, subset)).distinct()
        if self.journals_filter_flag == 'only':
            if self.journals_filter:
                subset = []
                for journal, val in self.journals_filter:
                    if val:
                        subset.append(Q(paper__journal__title__icontains=journal))
                if subset:
                    queryset = queryset.filter(reduce(operator.or_, subset)).distinct()
                else:
                    queryset = None
            else:
                queryset = None

        # authors
        if self.authors_filter_flag == 'all':
            if self.authors_filter:
                subset = []
                for pk, val in self.authors_filter:
                    if not val:
                        subset.append(Q(paper__authors__pk=pk))
                if subset:
                    queryset = queryset.exclude(reduce(operator.or_, subset)).distinct()
        if self.authors_filter_flag == 'only':
            if self.authors_filter:
                subset = []
                for pk, val in self.authors_filter:
                    if val:
                        subset.append(Q(paper__authors__pk=pk))
                if subset:
                    queryset = queryset.filter(reduce(operator.or_, subset)).distinct()
                else:
                    queryset = None
            else:
                queryset = None

        # like filter
        if self.like_flag:
            like_pks = UserTaste.objects\
                .filter(user=self.request.user,
                        is_liked=True,
                        scoring_method=self.request.user.settings.scoring_method)\
                .values('paper__pk')
            queryset = queryset.filter(paper_id__in=like_pks)

        # search query
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            queryset = queryset.filter(Q(paper__title__icontains=q) |
                                   Q(paper__abstract__icontains=q) |
                                   Q(paper__journal__title__icontains=q) |
                                   Q(paper__authors__last_name__icontains=q) |
                                   Q(paper__authors__first_name__icontains=q))\
                .distinct()

        # sort queryset
        if self.sorting_flag == 'relevant':
            # this is by default for feed
            pass
        elif self.sorting_flag == 'recent':
            # order by date
            queryset = queryset.order_by(Coalesce('paper__date_ep',
                                                  'paper__date_pp',
                                                  'paper__date_fs').desc())
        elif self.sorting_flag == 'trendy':
            # order by altmetric score
            queryset = queryset.order_by('-paper__altmetric__score')
        elif self.sorting_flag == 'nothing':
            # let's shuffle the results
            queryset = queryset.order_by("?")

        return queryset

library = UserLibraryView.as_view()


class UserLibraryTrashView(UserLibraryView):

    def get_queryset(self):
        query_set = UserLibPaper.objects\
            .filter(userlib=self.request.user.lib,
                    is_trashed=True)
        return self.filter_queryset(query_set)

    def get_context_data(self, **kwargs):
        context = super(UserLibraryTrashView, self).get_context_data(**kwargs)
        # library tab
        context['tab'] = 'trash'
        return context

library_trash = UserLibraryTrashView.as_view()


class UserLibraryLikesView(UserLibraryView):

    def get_queryset(self):

        query_set = UserTaste.objects\
            .filter(user=self.request.user,
                    is_liked=True)

        return self.filter_queryset(query_set)

    def get_context_data(self, **kwargs):
        context = super(UserLibraryLikesView, self).get_context_data(**kwargs)
        # library tab
        context['tab'] = 'likes'
        return context

library_likes = UserLibraryLikesView.as_view()


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
        self.request.user.settings.stream_time_lapse = form.cleaned_data['stream_time_lapse']
        self.request.user.settings.stream_narrowness = form.cleaned_data['stream_narrowness']
        # self.request.user.settings.stream_model = form.cleaned_data['stream_model']
        # self.request.user.settings.stream_method = form.cleaned_data['stream_method']
        self.request.user.settings.save()
        return super(UserStreamSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'stream_narrowness': self.request.user.settings.get_stream_narrowness_display(),
                'stream_time_lapse': self.request.user.settings.get_stream_time_lapse_display(),
                # 'stream_method': self.request.user.settings.stream_method,
                # 'stream_model': self.request.user.settings.stream_model,
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
        self.request.user.settings.trend_narrowness = form.cleaned_data['trend_narrowness']
        self.request.user.settings.trend_time_lapse = form.cleaned_data['trend_time_lapse']
        # self.request.user.settings.trend_model = form.cleaned_data['trend_model']
        # self.request.user.settings.trend_method = form.cleaned_data['trend_method']
        self.request.user.settings.save()
        return super(UserTrendSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'trend_narrowness': self.request.user.settings.get_trend_narrowness_display(),
                # 'trend_model': self.request.user.settings.trend_model.name,
                'trend_time_lapse': self.request.user.settings.get_trend_time_lapse_display(),
                # 'trend_method': self.request.user.settings.trend_method,
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
def update_library(request):
    if request.is_ajax():
        if request.user.lib.state == 'IDL':
            provider_name = request.user.social_auth.first().provider
            update_lib.delay(request.user.pk, provider_name)
            data = {}
            return JsonResponse(data)
        else:
            data = {'error': 'library is already syncing'}
            return JsonResponse(data, status=400)


@login_required
def like_call(request):
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
def tick_call(request):
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

        subject = 'An invitation to try PubStream'
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