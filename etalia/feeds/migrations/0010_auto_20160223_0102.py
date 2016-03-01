# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0009_trend_state'),
    ]

    operations = [
        migrations.AddField(
            model_name='streammatches',
            name='seen',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='trendmatches',
            name='seen',
            field=models.BooleanField(default=False),
        ),
    ]
