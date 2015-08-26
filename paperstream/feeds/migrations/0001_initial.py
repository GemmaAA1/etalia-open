# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserFeed',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=100)),
                ('state', models.CharField(blank=True, default='NON', max_length=3, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')])),
                ('message', models.CharField(default='', max_length=127, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserFeedMatchPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('is_score_computed', models.BooleanField(default=False)),
                ('is_disliked', models.BooleanField(default=False)),
                ('is_liked', models.BooleanField(default=False)),
                ('feed', models.ForeignKey(to='feeds.UserFeed')),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.CreateModel(
            name='UserFeedSeedPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('feed', models.ForeignKey(to='feeds.UserFeed')),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='UserFeedVector',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(null=True, base_field=models.FloatField(null=True), size=300)),
                ('feed', models.ForeignKey(related_name='vectors', to='feeds.UserFeed')),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
        ),
        migrations.AddField(
            model_name='userfeed',
            name='papers_match',
            field=models.ManyToManyField(related_name='feed_match', to='library.Paper', through='feeds.UserFeedMatchPaper'),
        ),
        migrations.AddField(
            model_name='userfeed',
            name='papers_seed',
            field=models.ManyToManyField(related_name='feed_seed', to='library.Paper', through='feeds.UserFeedSeedPaper'),
        ),
        migrations.AddField(
            model_name='userfeed',
            name='user',
            field=models.ForeignKey(related_name='feeds', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AlterUniqueTogether(
            name='userfeedvector',
            unique_together=set([('feed', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='userfeedseedpaper',
            unique_together=set([('feed', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='userfeedmatchpaper',
            unique_together=set([('feed', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='userfeed',
            unique_together=set([('name', 'user')]),
        ),
    ]
