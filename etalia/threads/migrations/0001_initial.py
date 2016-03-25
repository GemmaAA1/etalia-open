# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('library', '0003_journal_lib_size'),
    ]

    operations = [
        migrations.CreateModel(
            name='Thread',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('type', models.IntegerField(default=1, choices=[(1, 'Question'), (2, 'Paper')])),
                ('title', models.CharField(max_length=256)),
                ('content', models.TextField(default='', null=True, blank=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='threads')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadFeed',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=128)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('status', models.CharField(default='NON', choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], blank=True, max_length=3)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadFeedThread',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('new', models.BooleanField(default=True)),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('thread_feed', models.ForeignKey(to='threads.ThreadFeed')),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.CreateModel(
            name='ThreadMember',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('first_joined_at', models.DateTimeField(auto_now_add=True)),
                ('last_left_at', models.DateTimeField(default=None, null=True, blank=True)),
                ('num_comments', models.PositiveIntegerField(default=0)),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='ThreadNeighbor',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, verbose_name='Days from right now', choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')])),
                ('neighbors', django.contrib.postgres.fields.ArrayField(null=True, base_field=models.IntegerField(null=True), blank=True, size=10)),
                ('thread', models.ForeignKey(to='threads.Thread')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadPost',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('position', models.PositiveIntegerField(default=0)),
                ('content', models.TextField(default='', blank=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('thread', models.ForeignKey(to='threads.Thread', related_name='posts')),
            ],
        ),
        migrations.CreateModel(
            name='ThreadPostComment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('position', models.PositiveIntegerField(default=0)),
                ('content', models.TextField(default='', blank=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('post', models.ForeignKey(to='threads.ThreadPost', related_name='comments')),
            ],
        ),
        migrations.CreateModel(
            name='ThreadUserInvite',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('status', models.IntegerField(default=1, choices=[(1, 'Pending'), (2, 'Accepted'), (3, 'Declined')])),
                ('from_user', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='invite_from_users')),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('to_user', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='invite_to_users')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadUserState',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('is_pinned', models.BooleanField(default=False)),
                ('is_banned', models.BooleanField(default=False)),
                ('is_added', models.BooleanField(default=False)),
                ('is_trashed', models.BooleanField(default=False)),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='threads',
            field=models.ManyToManyField(to='threads.Thread', through='threads.ThreadFeedThread'),
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='thread',
            name='members',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL, through='threads.ThreadMember'),
        ),
        migrations.AddField(
            model_name='thread',
            name='paper',
            field=models.ForeignKey(to='library.Paper', null=True, blank=True, default=None),
        ),
        migrations.AlterUniqueTogether(
            name='threadpostcomment',
            unique_together=set([('post', 'position')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadpost',
            unique_together=set([('thread', 'position')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadfeedthread',
            unique_together=set([('thread', 'thread_feed')]),
        ),
    ]
