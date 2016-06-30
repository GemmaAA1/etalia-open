# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0010_paperengine_score_n_papers'),
    ]

    operations = [
        migrations.AddField(
            model_name='paperengine',
            name='score_n_threads',
            field=models.PositiveIntegerField(default=3650),
        ),
        migrations.AlterField(
            model_name='paperengine',
            name='score_n_papers',
            field=models.PositiveIntegerField(default=3650),
        ),
    ]
