from django.views.generic.base import ContextMixin

from .forms import CreateUserFeedForm


class CreateFeedModalMixin(ContextMixin):
    """Mixin to manage form in profile modal"""

    def get_context_data(self, **kwargs):
        context = super(CreateFeedModalMixin, self).get_context_data(**kwargs)
        # data = self.request.POST.copy()
        # data['user'] = self.request.user
        context['form_create_feed'] = CreateUserFeedForm()
        return context
