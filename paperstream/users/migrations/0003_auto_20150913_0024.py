# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0002_auto_20150909_1819'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='time_lapse',
            field=models.IntegerField(verbose_name='In the past for', choices=[(7, '1 Week'), (30, '1 Month'), (60, '2 Months'), (-1, 'All')], default=30),
        ),
    ]
