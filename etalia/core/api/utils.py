# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework.renderers import BrowsableAPIRenderer


class BrowsableAPIRendererWithoutForms(BrowsableAPIRenderer):
    """Renders the browsable api, but excludes the forms."""

    def get_context(self, *args, **kwargs):
        ctx = super(BrowsableAPIRendererWithoutForms, self).get_context(*args, **kwargs)
        ctx['display_edit_forms'] = False
        return ctx