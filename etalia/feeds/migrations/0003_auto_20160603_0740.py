# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0003_auto_20160603_0740'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('feeds', '0002_auto_20160602_0013'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadFeed',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(max_length=128, default='main')),
                ('state', models.CharField(blank=True, max_length=3, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')])),
                ('last_update', models.DateTimeField(blank=True, default=None, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='ThreadFeedThreads',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('date', models.DateField()),
                ('new', models.BooleanField(default=True)),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('threadfeed', models.ForeignKey(to='feeds.ThreadFeed')),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='threads',
            field=models.ManyToManyField(to='threads.Thread', through='feeds.ThreadFeedThreads'),
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='threadfeeds'),
        ),
        migrations.AlterUniqueTogether(
            name='threadfeedthreads',
            unique_together=set([('threadfeed', 'thread')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadfeed',
            unique_together=set([('user', 'name')]),
        ),
    ]
