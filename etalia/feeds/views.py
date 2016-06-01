# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.http import JsonResponse
from django.core.urlresolvers import reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings

from braces.views import LoginRequiredMixin

from etalia.users.models import UserTaste
from etalia.core.views import BasePaperListView
from etalia.core.mixins import NavFlapMixin, XMLMixin

from .models import Stream, StreamMatches, TrendMatches
from .tasks import update_stream, update_trend, reset_stream, reset_trend


class FeedPaperListView(BasePaperListView):

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.stream_time_lapse,
            'method': self.request.user.settings.stream_method,
        }
        return self.context_settings

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

        # filter time span
        self.original_qs = self.filter_time_span(self.original_qs)

        # filter
        queryset = self.original_qs
        if self.journals_filter:
            queryset = self.filter_journals(queryset)

        if self.authors_filter:
            queryset = self.filter_authors(queryset)

        if self.like_flag:
            queryset = self.filter_pin(queryset)

        if self.search_query:
            queryset = self.filter_search_query(queryset)

        return queryset

    def get_context_data(self, **kwargs):
        context = super(FeedPaperListView, self).get_context_data(**kwargs)
        context.update(self.get_context_endless(**kwargs))
        context.update(self.get_context_stats())
        context.update(self.get_context_usertaste())
        context.update(self.get_context_userlib())
        context.update(self.get_context_journals_filter())
        context.update(self.get_context_authors_filter())
        context.update(self.get_context_time_span())
        context.update(self.get_context_search_query())
        context.update(self.get_context_new_objects_since_last_login())
        context.update(self.get_context_control_session())

        return context

    def get_original_queryset(self):
        raise NotImplemented


class BaseStreamView(FeedPaperListView):
    model = StreamMatches
    template_name = 'feeds/feed.html'
    control_session = 'control_stream'

    def get_original_queryset(self):
        return self.model.objects\
            .filter(stream__name=self.kwargs.get('name', 'main'),
                    stream__user=self.request.user)\
            .select_related('paper',
                            'paper__journal',
                            'paper__altmetric')


class StreamView(LoginRequiredMixin, NavFlapMixin, BaseStreamView):
    page_template = 'feeds/feed_sub_page.html'

stream = StreamView.as_view()


class StreamViewXML(LoginRequiredMixin, NavFlapMixin, XMLMixin, BaseStreamView):
    page_template = 'feeds/feed_sub_page2.html'

stream_xml = StreamViewXML.as_view()


class BaseTrendView(FeedPaperListView):
    model = TrendMatches
    template_name = 'feeds/trend.html'
    control_session = 'control_trend'

    def get_original_queryset(self):
        return self.model.objects\
            .filter(trend__name=self.kwargs.get('name', 'main'),
                    trend__user=self.request.user)\
            .select_related('paper',
                            'paper__journal',
                            'paper__altmetric')


class TrendView(LoginRequiredMixin, NavFlapMixin, BaseTrendView):
    page_template = 'feeds/trend_sub_page.html'

trend = TrendView.as_view()


class TrendViewXML(LoginRequiredMixin, NavFlapMixin, XMLMixin, BaseTrendView):
    page_template = 'feeds/trend_sub_page2.html'

trend_xml = TrendViewXML.as_view()


@login_required
def reset_stream_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        stream_name = name
        reset_stream.delay(request.user.pk, stream_name=stream_name,
                           restrict_journal=False)
        data = {'display_update_modal': True,
                'message': 'Stream reset launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_stream_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        stream_name = name
        update_stream.delay(request.user.pk, stream_name=stream_name)
        data = {'display_update_modal': True,
                'message': 'Stream update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def update_trend_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        trend_name = name
        update_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'display_update_modal': True,
                'message': 'Trend update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def reset_trend_view(request, name):
    if request.is_ajax() or settings.DEBUG:
        trend_name = name
        reset_trend.delay(request.user.pk, trend_name=trend_name)
        data = {'display_update_modal': True,
                'message': 'Trend update launched.'}
        return JsonResponse(data)
    else:
        return redirect('feeds:stream')


@login_required
def ajax_user_feed_message(request, name):
    if request.method == 'GET':
        userfeed = get_object_or_404(Stream,
                                     name=name,
                                     user=request.user)
        if userfeed.state == 'IDL':
            data = {'done': True,
                    'url': str(reverse_lazy('feeds:feed',
                                        kwargs={'name': userfeed.name}))}
        else:
            data = {'done': False,
                    'message': userfeed.message}
        return JsonResponse(data)