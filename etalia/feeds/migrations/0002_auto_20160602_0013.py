# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='trend',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='trends'),
        ),
        migrations.AddField(
            model_name='streampapers',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AddField(
            model_name='streampapers',
            name='stream',
            field=models.ForeignKey(to='feeds.Stream'),
        ),
        migrations.AddField(
            model_name='stream',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='feeds.StreamPapers', related_name='stream_papers'),
        ),
        migrations.AddField(
            model_name='stream',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='streams'),
        ),
        migrations.AlterUniqueTogether(
            name='trendpapers',
            unique_together=set([('trend', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='trend',
            unique_together=set([('user', 'name')]),
        ),
        migrations.AlterUniqueTogether(
            name='streampapers',
            unique_together=set([('stream', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='stream',
            unique_together=set([('name', 'user')]),
        ),
    ]
