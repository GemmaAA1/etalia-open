from django.conf import settings
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.views.generic.list import ListView
from django.db.models import Q

from braces.views import LoginRequiredMixin

from core.mixins import ProfileModalFormsMixin
from library.models import Paper
from .models import UserFeed, UserFeedPaper


class home_feed(LoginRequiredMixin, ProfileModalFormsMixin, ListView):
    model = UserFeedPaper
    paginate_by = 10
    template_name = 'feeds/feed.html'

    def get_queryset(self):
        ufp = UserFeedPaper.objects.filter(
            Q(feed=self.request.user.feed.first()),
            Q(is_disliked=False))[:100]
        return ufp

home = home_feed.as_view()


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