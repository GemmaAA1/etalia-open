# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0006_auto_20151009_2004'),
    ]

    operations = [
        migrations.AddField(
            model_name='discoverfeed',
            name='time_lapse_top_altmetric',
            field=models.IntegerField(default=15),
        ),
        migrations.AlterField(
            model_name='discoverfeed',
            name='time_lapse',
            field=models.IntegerField(default=60),
        ),
    ]
