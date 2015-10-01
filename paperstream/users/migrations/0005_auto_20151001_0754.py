# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0004_auto_20150923_1734'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='scoring_method',
            field=models.IntegerField(verbose_name='Scoring Algo', default=1, choices=[(0, 'Max'), (1, 'Average'), (2, 'Average Threshold binary'), (3, 'Average weighted journal'), (4, 'Average weighted journal-date')]),
        ),
    ]
