# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0002_auto_20160602_0013'),
    ]

    operations = [
        migrations.AddField(
            model_name='paperengine',
            name='score_altmetric_boost',
            field=models.FloatField(default=0.1),
        ),
    ]
