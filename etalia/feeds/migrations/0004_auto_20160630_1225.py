# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0003_auto_20160603_0740'),
    ]

    operations = [
        migrations.AddField(
            model_name='stream',
            name='score_threshold',
            field=models.FloatField(default=0.6),
        ),
        migrations.AddField(
            model_name='trend',
            name='score_threshold',
            field=models.FloatField(default=0.6),
        ),
    ]
