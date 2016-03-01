# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0010_auto_20160223_0102'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='streammatches',
            name='seen',
        ),
        migrations.RemoveField(
            model_name='trendmatches',
            name='seen',
        ),
        migrations.AddField(
            model_name='streammatches',
            name='new',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='trendmatches',
            name='new',
            field=models.BooleanField(default=True),
        ),
    ]
