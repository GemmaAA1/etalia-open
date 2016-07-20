# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0005_threadfeed_score_threshold'),
    ]

    operations = [
        migrations.AlterField(
            model_name='stream',
            name='score_threshold',
            field=models.FloatField(default=0.3),
        ),
        migrations.AlterField(
            model_name='trend',
            name='score_threshold',
            field=models.FloatField(default=0.1),
        ),
    ]
