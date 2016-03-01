# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0008_auto_20160119_2328'),
    ]

    operations = [
        migrations.AddField(
            model_name='trend',
            name='state',
            field=models.CharField(choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], blank=True, max_length=3),
        ),
    ]
