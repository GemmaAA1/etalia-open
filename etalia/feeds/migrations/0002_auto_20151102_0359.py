# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0001_initial'),
        ('nlp', '0001_initial'),
        ('library', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.AddField(
            model_name='trend',
            name='user',
            field=models.ForeignKey(related_name='trend', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='streamvector',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='streamvector',
            name='stream',
            field=models.ForeignKey(related_name='vectors', to='feeds.Stream'),
        ),
        migrations.AddField(
            model_name='streamseeds',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AddField(
            model_name='streamseeds',
            name='stream',
            field=models.ForeignKey(to='feeds.Stream'),
        ),
        migrations.AddField(
            model_name='streammatches',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AddField(
            model_name='streammatches',
            name='stream',
            field=models.ForeignKey(to='feeds.Stream'),
        ),
        migrations.AddField(
            model_name='stream',
            name='matches',
            field=models.ManyToManyField(through='feeds.StreamMatches', related_name='stream_matches', to='library.Paper'),
        ),
        migrations.AddField(
            model_name='stream',
            name='seeds',
            field=models.ManyToManyField(through='feeds.StreamSeeds', related_name='stream_seeds', to='library.Paper'),
        ),
        migrations.AddField(
            model_name='stream',
            name='user',
            field=models.ForeignKey(related_name='streams', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AlterUniqueTogether(
            name='trendmatches',
            unique_together=set([('trend_feed', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='trend',
            unique_together=set([('user', 'name')]),
        ),
        migrations.AlterUniqueTogether(
            name='streamvector',
            unique_together=set([('stream', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='streamseeds',
            unique_together=set([('stream', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='streammatches',
            unique_together=set([('stream', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='stream',
            unique_together=set([('name', 'user')]),
        ),
    ]
