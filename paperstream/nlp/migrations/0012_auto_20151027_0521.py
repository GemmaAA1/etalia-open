# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0011_auto_20151021_0002'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperneighbors',
            name='time_lapse',
            field=models.IntegerField(verbose_name='Days from right now', choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (365, 'Year'), (-1, 'All')], default=-1),
        ),
    ]
