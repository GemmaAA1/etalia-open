# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0003_auto_20150521_2306'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='feed_status',
            field=models.CharField(max_length=3, default='', choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')]),
        ),
        migrations.AlterField(
            model_name='user',
            name='library_status',
            field=models.CharField(max_length=3, default='', choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')]),
        ),
    ]
