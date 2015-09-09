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
            name='time_lapse',
            field=models.IntegerField(verbose_name='In the past for', choices=[(7, '1 Week'), (30, '1 Month'), (60, '2 Months'), (-1, 'All')], default=61),
        ),
    ]
