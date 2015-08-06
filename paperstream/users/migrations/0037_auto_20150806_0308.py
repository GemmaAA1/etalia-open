# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0036_auto_20150805_0704'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='time_lapse',
            field=models.IntegerField(default=61, verbose_name='In the past for', choices=[(61, '2 Months'), (-1, 'All')]),
        ),
    ]
