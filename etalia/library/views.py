# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.template.response import TemplateResponse
from django.views.generic import DetailView, RedirectView
from django_filters.views import FilterMixin, FilterView
from django.utils.text import slugify

from .models import Paper
from .filters import PaperFilter, MyPaperFilter
from .mixins import PaperEagerLoadingMixin


def my_papers(request):
    return TemplateResponse(
        request,
        'papers/my_list.html',
        {'control_states': json.dumps(
            request.session.get('library-control-states',
                                {'time_span': None,
                                 'search': None,
                                 'pin': 0}
                                )
        )
        }
    )


class PaperDetail(DetailView):

    template_name = 'papers/detail.html'
    model = Paper

paper_slug = PaperDetail.as_view()


class PaperViewPk(RedirectView):
    """Redirect to slug paper url"""

    permanent = False
    query_string = True
    pattern_name = 'library:paper-slug'

    def get_redirect_url(self, *args, **kwargs):
        paper = Paper.objects.get(pk=kwargs['pk'])
        kwargs['slug'] = slugify(paper.title)
        return super(PaperViewPk, self).get_redirect_url(*args, **kwargs)

paper = PaperViewPk.as_view()


class PaperListView(PaperEagerLoadingMixin, FilterView):

    template_name = 'papers/list.html'
    # filterset_class = PaperFilter
    filterset_class = MyPaperFilter
    paginate_by = 25

    def get_filterset(self, filterset_class):
        filterset = super(PaperListView, self).get_filterset(filterset_class)
        filterset._qs = self.setup_eager_loading(filterset.qs,
                                                 user=self.request.user,
                                                 request=self.request)
        return filterset

papers = PaperListView.as_view()
