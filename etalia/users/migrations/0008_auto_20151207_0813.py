# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0007_auto_20151207_0607'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_method',
            field=models.IntegerField(default=0, choices=[(0, 'Occurrence Count'), (1, 'Max'), (2, 'Average'), (3, 'Average Threshold binary'), (4, 'Average weighted journal'), (5, 'Average weighted journal-date')], verbose_name='Method'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_method',
            field=models.IntegerField(default=0, choices=[(0, 'Method #1')], verbose_name='Method'),
        ),
    ]
