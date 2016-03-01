# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0008_auto_20151207_0813'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_narrowness',
            field=models.IntegerField(default=1, choices=[(-1, 'Narrower'), (0, 'Normal'), (1, 'Broader')], verbose_name='Narrowness'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_time_lapse',
            field=models.IntegerField(default=60, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Time range'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_narrowness',
            field=models.IntegerField(default=1, choices=[(-1, 'Narrower'), (0, 'Normal'), (1, 'Broader')], verbose_name='Narrowness'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_time_lapse',
            field=models.IntegerField(default=60, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Time range'),
        ),
    ]
