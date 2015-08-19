from django.views.generic.base import ContextMixin

from .forms import UpdateUserBasicForm, UserAffiliationForm, \
    UserSettingsForm


class ProfileModalFormsMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(ProfileModalFormsMixin, self).get_context_data(**kwargs)
        context['form_userbasic'] = \
            UpdateUserBasicForm(instance=self.request.user)
        context['form_affiliation'] = \
            UserAffiliationForm(instance=self.request.user.affiliation)
        context['form_settings'] = \
            UserSettingsForm(instance=self.request.user.settings)
        return context
