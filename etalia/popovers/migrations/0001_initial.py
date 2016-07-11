# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='PopOver',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('title', models.CharField(max_length=256)),
                ('body', models.TextField()),
                ('anchor', models.CharField(max_length=128)),
                ('type', models.PositiveIntegerField(default=1, choices=[(1, 'Anchored'), (2, 'Modal')])),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='UserPopOver',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('status', models.PositiveIntegerField(default=1, choices=[(1, 'New'), (2, 'Got It')])),
                ('popover', models.ForeignKey(to='popovers.PopOver')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
