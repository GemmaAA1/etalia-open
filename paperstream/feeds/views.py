from django.shortcuts import render
from braces.views import LoginRequiredMixin
from django.views.generic.list import ListView

from library.models import Paper

class home_feed(LoginRequiredMixin, ListView):
    model = Paper
    template_name = 'feeds/feed.html'

    def get_queryset(self):
        return Paper.objects.all()

home = home_feed.as_view()



@login_required
def dislikes(request):

    fnp_id = None
    if request.method == 'GET':
        fnp_id = request.GET['feednewpaper_id']

    if fnp_id:
        fnp = FeedNewPaper.objects.get(id=int(fnp_id))
        if fnp:
            fnp.is_disliked = not fnp.is_disliked
            if fnp.is_disliked:
                fnp.is_liked = False
            fnp.save()

        return HttpResponse(json.dumps({'is_disliked': fnp.is_disliked,
                                        'fnp_id': fnp_id}),
                            content_type="application/json")