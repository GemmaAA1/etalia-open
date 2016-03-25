# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.views.generic import UpdateView, FormView, DetailView, CreateView, \
    DeleteView, ListView
from braces.views import LoginRequiredMixin, UserPassesTestMixin
from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404

from etalia.core.mixins import AjaxableResponseMixin
from etalia.users.models import UserThread
from .forms import ThreadCreateForm, ThreadUpdateForm, ThreadPostForm, \
    ThreadPostCommentForm, UserThreadForm
from .models import Thread, ThreadPost, ThreadPostComment

from .serializers import ThreadSerializer, ThreadPostSerializer, \
    ThreadPostCommentSerializer, UserThreadSerializer


class ThreadView(LoginRequiredMixin, AjaxableResponseMixin, DetailView):
    model = Thread
    template_name = 'threads/thread.html'
    max_members = 3
    context_object_name = 'thread'

    def get_context_data(self, **kwargs):
        context = super(ThreadView, self).get_context_data(**kwargs)
        following = self.request.user.get_following()
        members = self.object.get_active_members()
        if following:
            members_non_followed = members.exclude(id=following)
        else:
            members_non_followed = members
        posts = self.object.posts.all()

        context['nb_members'] = members.count()
        context['nb_posts'] = posts.count()
        context['following'] = following
        context['state'], _ = UserThread.objects.get_or_create(
            user=self.request.user, thread=self.object)

        if context['state'].is_joined:
            context['members'] = members
            context['posts'] = posts
        else:
            context['members'] = []
            # restrict members to max_members, order by following + nb_comments
            context['members'] = [m for m in members if m in following][
                                 :self.max_members]
            context['members'][len(context['members']):self.max_members] = \
                members_non_followed[:self.max_members - len(context['members'])]

        return context

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': ThreadContextualizedSerializer(
                instance=self.object,
                context={'request': self.request}).data,
        }


thread = ThreadView.as_view()


class ThreadCreate(LoginRequiredMixin, AjaxableResponseMixin, CreateView):
    template_name = 'threads/thread_new.html'
    form_class = ThreadCreateForm

    def get_success_url(self, **kwargs):
        return reverse('threads:thread', kwargs={'pk': self.object.id})

    def get_form_kwargs(self):
        kwargs = super(ThreadCreate, self).get_form_kwargs()
        kwargs['user_id'] = self.request.user.id
        kwargs['paper_qs'] = self.request.user.lib.papers.all()
        return kwargs

    def form_valid(self, form):
        form.instance.owner_id = self.request.user.id
        return super(ThreadCreate, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        return {
            'redirect': self.get_success_url(),
            'results': ThreadSerializer(instance=self.object,
                                        context={'request': self.request}).data,
        }


new_thread = ThreadCreate.as_view()


class ThreadUpdateView(LoginRequiredMixin, UserPassesTestMixin,
                       AjaxableResponseMixin, UpdateView):
    form_class = ThreadUpdateForm
    template_name = 'threads/thread_update.html'

    def get_object(self, queryset=None):
        return get_object_or_404(Thread, pk=self.kwargs.get('pk'))

    def test_func(self, user):
        # test if user is thread owner
        try:
            if Thread.objects.get(pk=self.kwargs.get('pk')).owner == user:
                return True
            else:
                return False
        except Thread.DoesNotExist:
            return False

    def get_form_kwargs(self):
        kwargs = super(ThreadUpdateView, self).get_form_kwargs()
        kwargs['user_id'] = self.request.user.id
        return kwargs

    def get_context_data(self, **kwargs):
        context = super(ThreadUpdateView, self).get_context_data(**kwargs)
        context['thread'] = self.object
        return context

    def get_success_url(self, **kwargs):
        return reverse('threads:thread', kwargs={'pk': self.object.id})

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': ThreadSerializer(instance=self.object,
                                        context={'request': self.request}).data,
        }


update_thread = ThreadUpdateView.as_view()


class ThreadPostCreateView(LoginRequiredMixin, UserPassesTestMixin,
                           AjaxableResponseMixin,
                           CreateView):
    template_name = 'threads/post_new.html'
    form_class = ThreadPostForm
    thread = None

    def dispatch(self, request, *args, **kwargs):
        self.thread = get_object_or_404(Thread, pk=self.kwargs.get('pk'))
        return super(ThreadPostCreateView, self).dispatch(request, *args,
                                                          **kwargs)

    def test_func(self, user):
        # test if user is in thread members
        if user in self.thread.get_active_members():
            return True
        else:
            return False

    def form_valid(self, form):
        form.instance.author_id = self.request.user.id
        form.instance.thread_id = self.thread.id
        return super(ThreadPostCreateView, self).form_valid(form)

    def get_context_data(self, **kwargs):
        context = super(ThreadPostCreateView, self).get_context_data(**kwargs)
        context['thread'] = self.thread
        return context

    def get_success_url(self, **kwargs):
        return reverse('threads:thread', kwargs={'pk': self.thread.id})

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': ThreadPostSerializer(instance=self.object).data,
        }


new_post = ThreadPostCreateView.as_view()


class ThreadPostUpdateView(LoginRequiredMixin, UserPassesTestMixin,
                           AjaxableResponseMixin, UpdateView):
    template_name = 'threads/post_update.html'
    form_class = ThreadPostForm
    context_object_name = 'post'

    def test_func(self, user):
        # test if user is author of post
        try:
            tp = ThreadPost.objects.get(pk=self.kwargs.get('pk'))
            if user.id == tp.author.id:
                return True
            else:
                return False
        except ThreadPost.DoesNotExist:
            return False

    def get_success_url(self, **kwargs):
        return reverse('threads:thread', kwargs={'pk': self.object.thread.id})

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': ThreadPostSerializer(instance=self.object).data,
        }


edit_post = ThreadPostUpdateView.as_view()


class ThreadPostDeleteView(LoginRequiredMixin, AjaxableResponseMixin,
                           DeleteView):
    model = ThreadPost

    def get_success_url(self):
        return reverse('threads:thread', kwargs={'pk': self.object.thread.id})

    def get_ajax_data(self, *args, **kwargs):
        return {
        }


delete_post = ThreadPostDeleteView.as_view()


class ThreadPostCommentCreateView(LoginRequiredMixin, UserPassesTestMixin,
                                  AjaxableResponseMixin, CreateView):
    template_name = 'threads/comment_new.html'
    form_class = ThreadPostCommentForm
    post_ = None

    def dispatch(self, request, *args, **kwargs):
        self.post_ = get_object_or_404(ThreadPost, pk=self.kwargs.get('pk'))
        return super(ThreadPostCommentCreateView, self).dispatch(request, *args,
                                                                 **kwargs)

    def test_func(self, user):
        # test if user is in thread members
        if user in self.post_.thread.get_active_members():
            return True
        else:
            return False

    def get_success_url(self, **kwargs):
        return reverse('threads:thread',
                       kwargs={'pk': self.object.post.thread.id})

    def form_valid(self, form):
        form.instance.author_id = self.request.user.id
        form.instance.post_id = self.post_.id
        return super(ThreadPostCommentCreateView, self).form_valid(form)

    def get_context_data(self, **kwargs):
        context = super(ThreadPostCommentCreateView, self).get_context_data(
            **kwargs)
        context['post'] = self.post_
        return context

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': ThreadPostCommentSerializer(instance=self.object).data,
        }


new_comment = ThreadPostCommentCreateView.as_view()


class ThreadPostCommentUpdateView(LoginRequiredMixin, UserPassesTestMixin,
                                  AjaxableResponseMixin, UpdateView):
    template_name = 'threads/comment_update.html'
    form_class = ThreadPostCommentForm
    context_object_name = 'comment'

    def test_func(self, user):
        # test if user is author of comment
        try:
            tpc = ThreadPostComment.objects.get(pk=self.kwargs.get('pk'))
            if user.id == tpc.author.id:
                return True
            else:
                return False
        except ThreadPostComment.DoesNotExist:
            return False

    def get_success_url(self, **kwargs):
        return reverse('threads:thread',
                       kwargs={'pk': self.object.post.thread.id})

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': ThreadPostCommentSerializer(instance=self.object).data,
        }


edit_comment = ThreadPostCommentUpdateView.as_view()


class ThreadPostCommentDeleteView(LoginRequiredMixin, AjaxableResponseMixin,
                                  DeleteView):
    model = ThreadPostComment

    def get_success_url(self):
        return reverse('threads:thread',
                       kwargs={'pk': self.object.post.thread.id})

    def get_ajax_data(self, *args, **kwargs):
        return {
        }


delete_comment = ThreadPostCommentDeleteView.as_view()


class MyThreadsView(LoginRequiredMixin, AjaxableResponseMixin, ListView):
    pass

my_threads = MyThreadsView.as_view()


class UserThreadView(LoginRequiredMixin, AjaxableResponseMixin, FormView):

    form_class = UserThreadForm
    thread_id = None
    action = None
    object = None

    def get_form_kwargs(self):
        kwargs = super(UserThreadView, self).get_form_kwargs()
        self.thread_id = kwargs['data']['thread']
        # copy data (post data is immutable)
        data = kwargs['data'].copy()
        # pop action (action not not an attribute of ThreadUserState)
        self.action = data['action']
        data.pop('action')
        return kwargs

    def get_success_url(self):
        return reverse('threads:thread', kwargs={'pk': self.object.thread_id})

    def form_valid(self, form):
        ut, _ = UserThread.objects.get_or_create(user=self.request.user,
                                                 thread_id=self.thread_id)
        if self.action == 'pin':
            ut.pin()
        if self.action == 'ban':
            ut.ban()
        if self.action == 'join':
            ut.join()
        if self.action == 'leave':
            ut.leave()
        self.object = ut
        return super(UserThreadView, self).form_valid(form)

    def get_ajax_data(self, *args, **kwargs):
        return {
            'results': UserThreadSerializer(instance=self.object).data,
        }

thread_state = UserThreadView.as_view()