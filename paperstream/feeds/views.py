from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.views.generic.list import ListView
from django.views.generic import FormView, UpdateView
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings
from django.db.models import Q

from braces.views import LoginRequiredMixin

from core.mixins import ModalMixin, AjaxableResponseMixin
from .models import UserFeed, UserFeedPaper
from .forms import CreateUserFeedForm, UpdateUserFeedForm
from .tasks import update_feed as async_update_feed

class FeedView(LoginRequiredMixin, ModalMixin, ListView):
    """ClassView for displaying a UserFeed instance"""
    model = UserFeedPaper
    paginate_by = 10
    template_name = 'feeds/feed.html'

    def get_queryset(self):
        super(FeedView, self).get_queryset()
        feed_name = self.kwargs.get('feed_name', 'main') or 'main'
        self.userfeed = get_object_or_404(UserFeed, user=self.request.user,
                                          name=feed_name)
        query_set = UserFeedPaper.objects.filter(feed=self.userfeed,
                                                 is_disliked=False)

        query_set = self.filter_queryset(query_set)
        return query_set[:settings.FEEDS_DISPLAY_N_PAPERS]

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(paper__title__icontains=q) |
                                   Q(paper__abstract__icontains=q) |
                                   Q(paper__journal__title__icontains=q) |
                                   Q(paper__authors__last_name__icontains=q) |
                                   Q(paper__authors__first_name__icontains=q))
        # No q is specified so we return queryset
        return queryset

    def get_context_data(self, **kwargs):
        context = super(FeedView, self).get_context_data(**kwargs)
        context['userfeed'] = self.userfeed
        return context

feed_view = FeedView.as_view()


class CreateFeedView(LoginRequiredMixin, ModalMixin, FormView):
    """ClassView to create a new UserFeed instance"""
    model = UserFeed
    form_class = CreateUserFeedForm

    def get_success_url(self):
        return reverse('feeds:modify-feed',
                       kwargs={'feed_name': self.object.name})

    def form_valid(self, form):
        self.object = self.model(**form.cleaned_data)
        self.object.save()
        return super(CreateFeedView, self).form_valid(form)

    def form_invalid(self, form):
        return super(CreateFeedView, self).form_invalid(form)

    def get_form_kwargs(self):
        """Return kwarg to be passed to Form (used to get user from Form)"""
        kwargs = super(CreateFeedView, self).get_form_kwargs()
        kwargs['request'] = self.request
        return kwargs

create_feed_view = CreateFeedView.as_view()


class UpdateFeedView(LoginRequiredMixin, ModalMixin, UpdateView):
    """ClassView to update a UserFeed instance"""

    model = UserFeed
    form_class = UpdateUserFeedForm
    template_name = 'feeds/update_feed.html'

    def get_object(self, queryset=None):
        return UserFeed.objects.get(user=self.request.user,
                                    name=self.kwargs['feed_name'])

    def get_success_url(self):
        return reverse('feeds:feed', kwargs={'feed_name': self.object.name})


modify_feed_view = UpdateFeedView.as_view()

@login_required
def update_feed_view(request, pk):
    # check if user is owner of feed
    user = UserFeed.objects.get(pk=pk).user
    assert user == request.user
    async_update_feed.delay(pk)


# @login_required
# def dislikes(request):
#
#     fnp_id = None
#     if request.method == 'GET':
#         fnp_id = request.GET['feednewpaper_id']
#
#     if fnp_id:
#         fnp = FeedNewPaper.objects.get(id=int(fnp_id))
#         if fnp:
#             fnp.is_disliked = not fnp.is_disliked
#             if fnp.is_disliked:
#                 fnp.is_liked = False
#             fnp.save()
#
#         return HttpResponse(json.dumps({'is_disliked': fnp.is_disliked,
#                                         'fnp_id': fnp_id}),
#                             content_type="application/json")