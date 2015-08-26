from django.http import JsonResponse
from paperstream.users.mixins import ProfileModalFormsMixin
from paperstream.feeds.mixins import CreateFeedModalMixin

class AjaxableResponseMixin(object):
    """
    Mixin to add AJAX support to a form.
    Must be used with an object-based FormView (e.g. UpdateView)
    """

    def form_invalid(self, form):
        response = super(AjaxableResponseMixin, self).form_invalid(form)
        if self.request.is_ajax():
            return JsonResponse(form.errors, status=400)
        else:
            return response

    def form_valid(self, form):
        # We make sure to call the parent's form_valid() method because
        # it might to do some processing (in the case of CreateView, it will
        # call form.save() for example)
        response = super(AjaxableResponseMixin, self).form_valid(form)
        if self.request.is_ajax():
            data = self.get_ajax_data()
            return JsonResponse(data)
        else:
            return response

    def get_ajax_data(self):
        raise NotImplementedError


class ModalMixin(ProfileModalFormsMixin, CreateFeedModalMixin):
    """Pull Mixin in one"""
    pass
