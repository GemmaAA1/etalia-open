# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0015_auto_20151026_2341'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='time_lapse',
            field=models.IntegerField(verbose_name='In the past for', choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (365, 'Year'), (-1, 'All')], default=30),
        ),
    ]
