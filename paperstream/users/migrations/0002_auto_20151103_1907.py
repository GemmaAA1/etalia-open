# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_method',
            field=models.IntegerField(verbose_name='Method', default=1, choices=[(0, 'Max'), (1, 'Average'), (2, 'Average Threshold binary'), (3, 'Average weighted journal'), (4, 'Average weighted journal-date')]),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_time_lapse',
            field=models.IntegerField(verbose_name='Time range', default=30, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (365, 'Year'), (-1, 'All')]),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_method',
            field=models.IntegerField(verbose_name='Method', default=1, choices=[(0, 'Method #1')]),
        ),
    ]
