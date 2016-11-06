# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import json
from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.views.generic import DetailView, RedirectView
from django.utils.text import slugify

from .models import Paper


@login_required()
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

    template_name = 'feeds/my_list.html'
    model = Paper

    def get_template_names(self):
        if self.request.user.is_anonymous():
            self.template_name = 'trends/list.html'
        return super(PaperDetail, self).get_template_names()

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


class PaperRedirectViewPk(RedirectView):

    permanent = True
    query_string = True
    pattern_name = 'library:paper-slug'

    def get_redirect_url(self, *args, **kwargs):
        paper = Paper.objects.get(pk=kwargs['pk'])
        kwargs['slug'] = slugify(paper.title)
        return super(PaperRedirectViewPk, self).get_redirect_url(*args, **kwargs)

paper_redirect = PaperRedirectViewPk.as_view()


class PaperRedirectSlugViewPk(RedirectView):

    permanent = True
    query_string = True
    pattern_name = 'library:paper-slug'

    def get_redirect_url(self, *args, **kwargs):
        return super(PaperRedirectSlugViewPk, self).get_redirect_url(*args, **kwargs)

paper_redirect_slug = PaperRedirectSlugViewPk.as_view()


def papers(request):
    return TemplateResponse(
        request,
        'trends/list.html',
        {'control_states': json.dumps(
            request.session.get('library-control-states',
                                {'time_span': None,
                                 'search': None,
                                 'pin': 0}
                                )
        )
        }
    )