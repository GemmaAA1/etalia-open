# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0006_auto_20160608_1959'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfingerprint',
            name='state',
            field=models.CharField(max_length=3, blank=True, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')]),
        ),
    ]
