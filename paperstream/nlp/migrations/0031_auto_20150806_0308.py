# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0030_auto_20150805_2355'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lsh',
            name='time_lapse',
            field=models.IntegerField(default=-1, verbose_name='Days from right now', choices=[(61, '2 Months'), (-1, 'All')]),
        ),
    ]
