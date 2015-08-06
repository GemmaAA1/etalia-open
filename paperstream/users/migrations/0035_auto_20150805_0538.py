# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0034_auto_20150805_0047'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userlib',
            name='state',
            field=models.CharField(max_length=3, blank=True, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], default='NON'),
        ),
    ]
