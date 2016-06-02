# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('threads', '0001_initial'),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='threaduserinvite',
            name='from_user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='invite_from_users'),
        ),
        migrations.AddField(
            model_name='threaduserinvite',
            name='thread',
            field=models.ForeignKey(to='threads.Thread'),
        ),
        migrations.AddField(
            model_name='threaduserinvite',
            name='to_user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='invite_to_users'),
        ),
        migrations.AddField(
            model_name='threaduserhistory',
            name='threaduser',
            field=models.ForeignKey(to='threads.ThreadUser'),
        ),
        migrations.AddField(
            model_name='threaduser',
            name='thread',
            field=models.ForeignKey(to='threads.Thread'),
        ),
        migrations.AddField(
            model_name='threaduser',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='threadscore',
            name='thread',
            field=models.ForeignKey(to='threads.Thread'),
        ),
        migrations.AddField(
            model_name='threadscore',
            name='thread_feed',
            field=models.ForeignKey(to='threads.ThreadFeed'),
        ),
        migrations.AddField(
            model_name='threadpost',
            name='thread',
            field=models.ForeignKey(to='threads.Thread', related_name='posts'),
        ),
        migrations.AddField(
            model_name='threadpost',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='thread_scores',
            field=models.ManyToManyField(to='threads.Thread', through='threads.ThreadScore'),
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='threadcomment',
            name='post',
            field=models.ForeignKey(to='threads.ThreadPost', related_name='comments'),
        ),
        migrations.AddField(
            model_name='threadcomment',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='thread',
            name='paper',
            field=models.ForeignKey(blank=True, default=None, to='library.Paper', null=True, verbose_name='Related Paper'),
        ),
        migrations.AddField(
            model_name='thread',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='threads_owned', null=True),
        ),
        migrations.AlterUniqueTogether(
            name='threaduserinvite',
            unique_together=set([('thread', 'from_user', 'to_user')]),
        ),
        migrations.AlterUniqueTogether(
            name='threaduser',
            unique_together=set([('thread', 'user')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadscore',
            unique_together=set([('thread', 'thread_feed')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadpost',
            unique_together=set([('thread', 'content', 'user')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadcomment',
            unique_together=set([('post', 'user', 'content')]),
        ),
    ]
