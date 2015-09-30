# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20150929_0710'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('feeds', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserTaste',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('is_disliked', models.BooleanField(default=False)),
                ('is_liked', models.BooleanField(default=False)),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='tastes')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='userfeedmatchpaper',
            name='is_disliked',
        ),
        migrations.RemoveField(
            model_name='userfeedmatchpaper',
            name='is_liked',
        ),
    ]
