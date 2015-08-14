# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0012_auto_20150813_1854'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userfeed',
            name='status',
        ),
        migrations.AddField(
            model_name='userfeed',
            name='state',
            field=models.CharField(blank=True, max_length=3, default='NON', choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')]),
        ),
    ]
