# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0021_auto_20150616_0508'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userstats',
            name='number_papers',
        ),
        migrations.AddField(
            model_name='userstats',
            name='options',
            field=models.CharField(default='', blank=True, max_length=64),
        ),
        migrations.AlterField(
            model_name='userstats',
            name='state',
            field=models.CharField(max_length=3, choices=[('LIN', 'Log in'), ('LOU', 'Log out'), ('LSS', 'Library starts syncing'), ('LES', 'Library ends syncing'), ('FSS', 'Feed starts sync'), ('FES', 'Feed ends sync'), ('EMA', 'Email validated'), ('CRE', 'Create user')]),
        ),
    ]
