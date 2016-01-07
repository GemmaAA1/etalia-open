# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0005_stream_last_update'),
    ]

    operations = [
        migrations.AddField(
            model_name='streammatches',
            name='date',
            field=models.DateField(null=True, blank=True),
        ),
    ]
