# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('library', '0007_auto_20151009_1916'),
        ('feeds', '0004_auto_20151008_0924'),
    ]

    operations = [
        migrations.CreateModel(
            name='DiscoverFeed',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=365)),
                ('size', models.IntegerField(default=400)),
                ('past_n_papers', models.IntegerField(default=200000)),
                ('top_n_closest', models.IntegerField(default=5000)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='DiscoverFeedPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('discover_feed', models.ForeignKey(to='feeds.DiscoverFeed')),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
        ),
        migrations.AddField(
            model_name='discoverfeed',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='feeds.DiscoverFeedPaper'),
        ),
        migrations.AddField(
            model_name='discoverfeed',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='discover'),
        ),
        migrations.AlterUniqueTogether(
            name='discoverfeedpaper',
            unique_together=set([('discover_feed', 'paper')]),
        ),
    ]
