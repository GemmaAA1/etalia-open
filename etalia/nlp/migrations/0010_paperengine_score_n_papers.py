# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0009_auto_20160627_1054'),
    ]

    operations = [
        migrations.AddField(
            model_name='paperengine',
            name='score_n_papers',
            field=models.PositiveIntegerField(default=1095),
        ),
    ]
