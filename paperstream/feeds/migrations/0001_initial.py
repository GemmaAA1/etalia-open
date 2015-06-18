# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_paper_id_isbn'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='UserFeed',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('status', models.CharField(default='', choices=[('', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], blank=True, max_length=3)),
            ],
        ),
        migrations.CreateModel(
            name='UserFeedPaper',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, primary_key=True, verbose_name='ID')),
                ('score', models.FloatField(default=0.0)),
                ('feed', models.ForeignKey(to='feeds.UserFeed')),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
        ),
        migrations.AddField(
            model_name='userfeed',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='feeds.UserFeedPaper'),
        ),
        migrations.AddField(
            model_name='userfeed',
            name='user',
            field=models.ForeignKey(related_name='feed', to=settings.AUTH_USER_MODEL),
        ),
    ]
