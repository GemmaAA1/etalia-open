# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0007_userfingerprint_state'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperengine',
            name='score_altmetric_boost',
            field=models.FloatField(default=0.4),
        ),
    ]
