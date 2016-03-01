# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0004_auto_20151105_0422'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperneighbors',
            name='time_lapse',
            field=models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Days from right now'),
        ),
    ]
