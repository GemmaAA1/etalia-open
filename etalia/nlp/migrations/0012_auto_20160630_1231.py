# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0011_auto_20160630_1230'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='paperengine',
            name='score_n_threads',
        ),
        migrations.AddField(
            model_name='threadengine',
            name='score_n_threads',
            field=models.PositiveIntegerField(default=3650),
        ),
    ]
