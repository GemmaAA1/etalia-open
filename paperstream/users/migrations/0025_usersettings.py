# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0010_auto_20150709_1953'),
        ('users', '0024_auto_20150708_0608'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserSettings',
            fields=[
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL, primary_key=True, serialize=False, related_name='settings')),
                ('time_lapse', models.IntegerField(choices=[(7, '1 Week'), (14, '2 Weeks'), (30, '1 Month'), (90, '3 Months'), (180, '6 Months'), (365, '1 Year'), (1825, '5 Year'), (-1, 'All')], default=7)),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
