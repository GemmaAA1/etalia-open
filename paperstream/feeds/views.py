# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.http import JsonResponse
from django.core.urlresolvers import reverse_lazy
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings

from paperstream.users.models import UserTaste
from paperstream.core.views import BasePaperListView

from .models import Stream, StreamMatches, TrendMatches
from .tasks import update_stream, update_trend, reset_stream, reset_trend


class BaseStreamView(BasePaperListView):
    model = StreamMatches
    template_name = 'feeds/feed.html'

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.stream_time_lapse,
            'method': self.request.user.settings.stream_method,
            'model': self.request.user.settings.stream_model,
        }
        return self.context_settings

    def get_queryset(self):

        # Retrieve get arguments if any
        self.parse_ajax_data()

        # Get data
        # get paper rejected
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')
        # get paper
        self.original_qs = self.model.objects\
            .filter(stream__name=self.kwargs.get('name', 'main'),
                    stream__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper',
                            'paper__journal')

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


class StreamView(BaseStreamView):
    page_template = 'feeds/feed_sub_page.html'

    def parse_ajax_data(self):
        if self.request.GET.dict().get('querystring_key'):  # endless scrolling
            self.return_filter = False
        super(StreamView, self).parse_ajax_data()

stream_view = StreamView.as_view()


class StreamViewXML(BaseStreamView):
    page_template = 'feeds/feed_sub_page2.html'

stream_view_xml = StreamViewXML.as_view()


class BaseTrendView(BasePaperListView):
    model = TrendMatches
    template_name = 'feeds/feed.html'

    def get_context_settings(self):
        self.context_settings = {
            'time_lapse': self.request.user.settings.trend_time_lapse,
            'method': self.request.user.settings.trend_method,
            'model': self.request.user.settings.trend_model,
        }
        return self.context_settings

    def get_queryset(self):

        self.parse_ajax_data()

        # get ticked paper
        papers_ticked = UserTaste.objects\
            .filter(user=self.request.user,
                    is_ticked=True)\
            .values('paper')

        self.original_qs = self.model.objects\
            .filter(trend__name=self.kwargs.get('name', 'main'),
                    trend__user=self.request.user)\
            .exclude(paper__in=papers_ticked)\
            .select_related('paper',
                            'paper__journal')

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


class TrendView(BaseTrendView):
    page_template = 'feeds/feed_sub_page.html'

trend_view = TrendView.as_view()


class TrendViewFilter(BaseTrendView):
    page_template = 'feeds/feed_sub_page2.html'

trend_view_xml = TrendViewFilter.as_view()


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