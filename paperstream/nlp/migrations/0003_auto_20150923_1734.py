# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0002_auto_20150909_1819'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lsh',
            name='time_lapse',
            field=models.IntegerField(verbose_name='Days from right now', default=-1, choices=[(7, '1 Week'), (30, '1 Month'), (60, '2 Months'), (365, '1 Year'), (-1, 'All')]),
        ),
    ]
