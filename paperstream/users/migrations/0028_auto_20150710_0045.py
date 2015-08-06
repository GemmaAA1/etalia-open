# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0027_auto_20150709_2319'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='model',
            field=models.ForeignKey(verbose_name='NLP Model', to='nlp.Model'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='time_lapse',
            field=models.IntegerField(verbose_name='Feed from the past', default=7, choices=[(7, '1 Week'), (14, '2 Weeks'), (30, '1 Month'), (90, '3 Months'), (180, '6 Months'), (365, '1 Year'), (1825, '5 Year'), (-1, 'All')]),
        ),
    ]
