# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('library', '0007_auto_20151009_1916'),
        ('feeds', '0007_auto_20151009_2331'),
    ]

    operations = [
        migrations.CreateModel(
            name='TrendFeed',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=60)),
                ('size', models.IntegerField(default=400)),
                ('past_n_papers', models.IntegerField(default=200000)),
                ('top_n_closest', models.IntegerField(default=5000)),
                ('time_lapse_top_altmetric', models.IntegerField(default=15)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='TrendFeedPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('trend_feed', models.ForeignKey(to='feeds.TrendFeed')),
            ],
        ),
        migrations.RemoveField(
            model_name='discoverfeed',
            name='matches',
        ),
        migrations.RemoveField(
            model_name='discoverfeed',
            name='user',
        ),
        migrations.AlterUniqueTogether(
            name='discoverfeedpaper',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='discoverfeedpaper',
            name='discover_feed',
        ),
        migrations.RemoveField(
            model_name='discoverfeedpaper',
            name='paper',
        ),
        migrations.DeleteModel(
            name='DiscoverFeed',
        ),
        migrations.DeleteModel(
            name='DiscoverFeedPaper',
        ),
        migrations.AddField(
            model_name='trendfeed',
            name='matches',
            field=models.ManyToManyField(to='library.Paper', through='feeds.TrendFeedPaper'),
        ),
        migrations.AddField(
            model_name='trendfeed',
            name='user',
            field=models.OneToOneField(to=settings.AUTH_USER_MODEL, related_name='trend'),
        ),
        migrations.AlterUniqueTogether(
            name='trendfeedpaper',
            unique_together=set([('trend_feed', 'paper')]),
        ),
    ]
