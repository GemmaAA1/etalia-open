from django.contrib.auth.decorators import login_required
from django.views.generic.list import ListView
from django.shortcuts import get_object_or_404
from django.conf import settings

from braces.views import LoginRequiredMixin

from core.mixins import ProfileModalFormsMixin
from .models import UserFeed, UserFeedPaper

from .tasks import update_feed as async_update_feed

class Feed(LoginRequiredMixin, ProfileModalFormsMixin, ListView):

    model = UserFeedPaper
    paginate_by = 10
    template_name = 'feeds/feed.html'

    def get_queryset(self):
        super(Feed, self).get_queryset()
        feed_name = self.kwargs.get('feed_name', 'main') or 'main'
        self.userfeed = get_object_or_404(UserFeed, user=self.request.user,
                                          name=feed_name)
        ufp = UserFeedPaper.objects.filter(feed=self.userfeed,
                                           is_disliked=False)
        return ufp[:settings.FEEDS_DISPLAY_N_PAPERS]

    def get_context_data(self, **kwargs):
        context = super(Feed, self).get_context_data(**kwargs)
        context['userfeed'] = self.userfeed
        return context

feed_view = Feed.as_view()

@login_required
def update_feed(request, pk):
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