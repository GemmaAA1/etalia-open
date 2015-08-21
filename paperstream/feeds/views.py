from django.http import JsonResponse
from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.views.generic.list import ListView
from django.views.generic import FormView, DeleteView, RedirectView
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings
from django.db.models import Q

from braces.views import LoginRequiredMixin

from core.mixins import ModalMixin, AjaxableResponseMixin
from library.models import Paper

from .models import UserFeed, UserFeedMatchPaper, UserFeedSeedPaper
from .forms import CreateUserFeedForm, CheckedForm
from .tasks import update_feed as async_update_feed


class FeedView(LoginRequiredMixin, ModalMixin, ListView):
    """ClassView for displaying a UserFeed instance"""
    model = Paper
    paginate_by = 10
    template_name = 'feeds/feed.html'

    def get_queryset(self):
        query_set = UserFeedMatchPaper.objects\
            .filter(feed=self.userfeed, is_disliked=False)\
            .select_related('paper')
        query_set = self.filter_queryset(query_set)
        return query_set

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(title__icontains=q) |
                                   Q(abstract__icontains=q) |
                                   Q(journal__title__icontains=q) |
                                   Q(authors__last_name__icontains=q) |
                                   Q(authors__first_name__icontains=q))\
                .distinct()
        # No q is specified so we return queryset
        return queryset

    def get_context_data(self, **kwargs):
        context = super(FeedView, self).get_context_data(**kwargs)
        context['userfeed'] = self.userfeed
        return context

    def get_userfeed(self, request, **kwargs):
        feed_pk = int(kwargs.get('pk', 0))
        if feed_pk:
            userfeed = get_object_or_404(UserFeed,
                                         pk=feed_pk,
                                         user=request.user)
        else:
            userfeed = get_object_or_404(UserFeed,
                                         user=request.user,
                                         name='main')
        return userfeed

    def get(self, request, *args, **kwargs):
        self.userfeed = self.get_userfeed(request, **kwargs)
        return super(FeedView, self).get(request, *args, **kwargs)

feed_view = FeedView.as_view()


class FeedMainView(LoginRequiredMixin, RedirectView):
    """Redirect to main feed"""

    pattern_name = 'feeds:feed'
    permanent = False

    def get_redirect_url(self, *args, **kwargs):
        # Check if userfeed pk matched db and current logged user
        userfeed = get_object_or_404(UserFeed, name='main',
                                     user=self.request.user)
        pk = userfeed.id
        return super(FeedMainView, self).get_redirect_url(pk=pk)

feed_main = FeedMainView.as_view()


class CreateFeedView(LoginRequiredMixin, ModalMixin, FormView):
    """ClassView to create a new UserFeed instance"""
    model = UserFeed
    form_class = CreateUserFeedForm

    def get_success_url(self):
        return reverse('feeds:modify-feed',
                       kwargs={'pk': self.object.pk})

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


class ModifyFeedView(LoginRequiredMixin, ModalMixin, ListView):

    model = Paper
    template_name = 'feeds/feed_settings.html'
    paginate_by = 50

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(title__icontains=q) |
                                   Q(abstract__icontains=q) |
                                   Q(journal__title__icontains=q) |
                                   Q(authors__last_name__icontains=q) |
                                   Q(authors__first_name__icontains=q))\
                .distinct()
        # No q is specified so we return queryset
        return queryset

    def get_queryset(self):
        # get all user lib paper
        queryset = self.request.user.lib.papers.all()
        # filter through search input
        queryset = self.filter_queryset(queryset)
        return queryset

    def get_context_data(self, **kwargs):
        context = super(ModifyFeedView, self).get_context_data(**kwargs)
        context['userfeed'] = self.userfeed
        context['seed_pks'] = self.userfeed.papers_seed.all()\
            .values_list('pk', flat=True)
        return context

    def get_userfeed(self, request, **kwargs):
        feed_pk = int(kwargs.get('pk', 0))
        if feed_pk:
            userfeed = get_object_or_404(UserFeed,
                                         pk=feed_pk,
                                         user=request.user)
        else:
            userfeed = get_object_or_404(UserFeed,
                                         user=request.user,
                                         name='main')
        return userfeed

    def get(self, request, *args, **kwargs):
        self.userfeed = self.get_userfeed(request, **kwargs)
        return super(ModifyFeedView, self).get(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        self.userfeed = self.get_userfeed(request, **kwargs)
        seed_pks = self.userfeed.papers_seed.all().values_list('pk', flat=True)
        if request.POST:
            pk = int(request.POST.get('pk'))
            if pk in seed_pks:      # remove it
                UserFeedSeedPaper.objects\
                    .get(feed=self.userfeed, paper_id=pk)\
                    .delete()
                data = {str(pk): ''}
            else:       # add it
                UserFeedSeedPaper.objects\
                    .create(feed=self.userfeed, paper_id=pk)
                data = {str(pk): 'success'}
            if request.is_ajax():
                return JsonResponse(data)
            else:
                return super(ModifyFeedView, self).get(self, request, *args,
                                                       **kwargs)

modify_feed_view = ModifyFeedView.as_view()


class DeleteFeedView(LoginRequiredMixin, ModalMixin, DeleteView):
    """ClassView to update a UserFeed instance"""

    model = UserFeed
    success_url = reverse_lazy('feeds:main')

delete_feed_view = DeleteFeedView.as_view()


class UpdateFeedView(LoginRequiredMixin, ModalMixin, RedirectView):

    pattern_name = 'feeds:feed'
    permanent = False

    def get_redirect_url(self, *args, **kwargs):
        # Check if userfeed pk matched db and current logged user
        userfeed = get_object_or_404(UserFeed,
                              pk=kwargs['pk'],
                              user=self.request.user)
        # update feed async
        userfeed.set_state('ING')
        async_update_feed.delay(kwargs['pk'])
        return super(UpdateFeedView, self).get_redirect_url(*args, **kwargs)

update_feed_view = UpdateFeedView.as_view()


@login_required
def ajax_user_feed_message(request, pk):
    if request.method == 'GET':
        userfeed = get_object_or_404(UserFeed,
                                     pk=pk,
                                     user=request.user)
        if userfeed.state == 'IDL':
            data = {'done': True,
                    'url': str(reverse_lazy('feeds:feed',
                                        kwargs={'pk': userfeed.id}))}
        else:
            data = {'done': False,
                    'message': userfeed.message}
        return JsonResponse(data)

@login_required
def feed_like_view(request, pk):
    if request.method == 'POST':
        ufp_pk = int(request.POST.get('pk'))
        ufmp = get_object_or_404(UserFeedMatchPaper, pk=ufp_pk)
        assert ufmp.feed.user == request.user
        if ufmp.is_liked:
            ufmp.is_liked = False
        else:
            ufmp.is_liked = True
            ufmp.is_disliked = False
        ufmp.save()
        data = {'is_liked': ufmp.is_liked,
                'is_disliked': ufmp.is_disliked}
        return JsonResponse(data)
    else:
        redirect('feeds:feed', kwargs={'pk': pk})


@login_required
def feed_dislike_view(request, pk):
    if request.method == 'POST':
        ufp_pk = int(request.POST.get('pk'))
        ufmp = get_object_or_404(UserFeedMatchPaper, pk=ufp_pk)
        assert ufmp.feed.user == request.user
        if ufmp.is_disliked:
            ufmp.is_disliked = False
        else:
            ufmp.is_liked = False
            ufmp.is_disliked = True
        ufmp.save()
        data = {'is_liked': ufmp.is_liked,
                'is_disliked': ufmp.is_disliked}
        return JsonResponse(data)
    else:
        redirect('feeds:feed', kwargs={'pk': pk})

