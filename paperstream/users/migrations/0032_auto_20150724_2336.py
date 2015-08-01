# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0031_usersettings_method'),
    ]

    operations = [
        migrations.RenameField(
            model_name='usersettings',
            old_name='method',
            new_name='scoring_method',
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='time_lapse',
            field=models.IntegerField(verbose_name='Feed from the past', default=7, choices=[(7, '1 Week'), (14, '2 Weeks'), (31, '1 Month'), (182, '6 Months'), (365, '1 Year')]),
        ),
    ]
