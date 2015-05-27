# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0004_auto_20150521_2307'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='feed_status',
            field=models.CharField(choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')], max_length=3, default='', blank=True),
        ),
        migrations.AlterField(
            model_name='user',
            name='library_status',
            field=models.CharField(choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')], max_length=3, default='', blank=True),
        ),
    ]
