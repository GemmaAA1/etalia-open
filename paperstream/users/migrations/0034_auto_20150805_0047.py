# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0033_auto_20150729_1817'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userlib',
            name='status',
        ),
        migrations.AddField(
            model_name='userlib',
            name='state',
            field=models.CharField(max_length=3, default='', blank=True, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')]),
        ),
    ]
