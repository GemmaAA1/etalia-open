# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0014_auto_20150528_0628'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserStats',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('state', models.CharField(max_length=3, choices=[('LIN', 'Log in'), ('LOU', 'Log out'), ('LIB', 'Library sync'), ('FEE', 'Feed sync')])),
                ('number_papers', models.IntegerField(default=0)),
                ('datetime', models.DateTimeField(auto_now_add=True)),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='stats')),
            ],
        ),
    ]
