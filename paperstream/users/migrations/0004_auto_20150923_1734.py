# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0003_auto_20150913_0024'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='time_lapse',
            field=models.IntegerField(verbose_name='In the past for', default=30, choices=[(7, '1 Week'), (30, '1 Month'), (60, '2 Months'), (365, '1 Year'), (-1, 'All')]),
        ),
    ]
