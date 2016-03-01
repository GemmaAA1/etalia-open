# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.http import Http404
from django.test import RequestFactory, Client
from django.conf import settings
from django.contrib.auth.models import AnonymousUser
from .base import UserFeedTestCase

from ..models import Stream
from ..views import feed_view, update_feed_view, create_feed_view

class FeedViewTestCase(UserFeedTestCase):

    def setUp(self):
        super(FeedViewTestCase, self).setUp()
        self.factory = RequestFactory()
        self.feed_default = Stream.objects.create_main(user=self.user)
        self.feed2 = Stream.objects.create(user=self.user, name='test',
                                             papers_seed=self.papers)

    def test_feed_page_redirect_to_login_if_user_anonymous(self):
        request = self.factory.get('/feed/test')
        request.user = AnonymousUser()
        response = feed_view(request, feed_name='main')
        self.assertEqual(response.status_code, 302)
        self.assertTrue(response.url.startswith(settings.LOGIN_REDIRECT_URL))

    def test_feed_page_returns_404_if_feed_name_is_unknown(self):
        pk_fake = 23143124123413241324
        request = self.factory.get('/feed/{pk}'
                                   .format(pk=pk_fake))
        request.user = self.user
        with self.assertRaises(Http404):
            feed_view(request, pk=pk_fake)

    def test_default_feed_page_renders_feed_template(self):
        request = self.factory.get('/feed')
        request.user = self.user
        response = feed_view(request)
        self.assertIn('feeds/feed.html', response.template_name)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context_data['userfeed'], self.feed_default)

    def test_feed_page_return_200(self):
        request = self.factory.get('/feed/{pk}'
                                   .format(pk=self.feed2.pk))
        request.user = self.user
        response = feed_view(request)
        self.assertEqual(response.status_code, 200)

    def test_feed_page_return(self):
        request = self.factory.get('/feed/{pk}'
                                   .format(pk=self.feed2.pk))
        request.user = self.user
        response = feed_view(request)
        self.assertEqual(response.status_code, 200)

    def test_feed_page_return_correct_userfeed(self):
        request = self.factory.get('/feed/{pk}'
                                   .format(pk=self.feed2.pk))
        request.user = self.user
        response = feed_view(request, pk=self.feed2.pk)
        self.assertEqual(response.context_data['userfeed'], self.feed2)

    def test_feed_page_list_object(self):
        request = self.factory.get('/feed/{pk}'
                                   .format(pk=self.feed2.pk))
        request.user = self.user
        response = feed_view(request, pk=self.feed2.pk)
        self.assertEqual(response.context_data['object_list'].count(),
                         min(self.feed2.papers_match.count(),
                             settings.FEEDS_DISPLAY_N_PAPERS))

class FeedUpdateTest(UserFeedTestCase):

    def setUp(self):
        super(FeedUpdateTest, self).setUp()
        self.factory = RequestFactory()
        self.feed_default = Stream.objects.create_main(user=self.user)
        self.feed2 = Stream.objects.create(user=self.user, name='test',
                                             papers_seed=self.papers)

    def test_user_can_run_update_feed(self):
        feed_pk = self.feed2.pk
        request = self.factory.get('/feed/{feed_pk}/update'
                                   .format(feed_pk=feed_pk))
        request.user = self.user
        update_feed_view(request, pk=feed_pk)


class CreateFeedTest(UserFeedTestCase):

    def setUp(self):
        super(CreateFeedTest, self).setUp()
        self.factory = RequestFactory()

    def test_create_feed_create_new_feed(self):
        feed_name = 'test'
        request = self.factory.post('/feed/create-feed',
                                    data={'name': feed_name})
        request.user = self.user
        create_feed_view(request)
        feed_names = Stream.objects.filter(user=self.user).values_list('name', flat=True)
        self.assertIn(feed_name, feed_names)

    def test_create_feed_cannot_create_feed_that_already_have_same_name(self):
        Stream.objects.create(user=self.user, name='test',
                                papers_seed=self.papers)
        feed_name = 'test'
        request = self.factory.post('/feed/create-feed',
                                    data={'name': feed_name})
        request.user = self.user
        response = create_feed_view(request)
        print('')