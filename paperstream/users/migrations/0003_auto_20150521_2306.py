# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0002_auto_20150521_2302'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='feed_status',
            field=models.CharField(max_length=3, choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')]),
        ),
        migrations.AlterField(
            model_name='user',
            name='library_status',
            field=models.CharField(max_length=3, choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')]),
        ),
    ]
