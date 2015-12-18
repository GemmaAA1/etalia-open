# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0004_auto_20151103_0109'),
    ]

    operations = [
        migrations.AddField(
            model_name='stream',
            name='last_update',
            field=models.DateTimeField(null=True, default=None, blank=True),
        ),
    ]
