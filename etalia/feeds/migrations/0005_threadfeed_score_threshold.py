# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0004_auto_20160630_1225'),
    ]

    operations = [
        migrations.AddField(
            model_name='threadfeed',
            name='score_threshold',
            field=models.FloatField(default=0.3),
        ),
    ]
