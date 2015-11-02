# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic.base import ContextMixin

from .forms import UpdateUserNameForm, UserAffiliationForm, \
    UserSettingsForm, UpdateUserTitleForm, UpdateUserPositionForm


class ProfileModalFormsMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(ProfileModalFormsMixin, self).get_context_data(**kwargs)
        if not self.request.user.is_anonymous():
            context['form_name'] = \
                UpdateUserNameForm(instance=self.request.user)
            context['form_title'] = \
                UpdateUserTitleForm(instance=self.request.user)
            context['form_position'] = \
                UpdateUserPositionForm(instance=self.request.user)
            context['form_affiliation'] = \
                UserAffiliationForm(instance=self.request.user.affiliation)
        return context


class SettingsModalFormsMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(SettingsModalFormsMixin, self).get_context_data(**kwargs)
        if not self.request.user.is_anonymous():
            context['form_settings'] = \
                UserSettingsForm(instance=self.request.user.settings)
        return context
