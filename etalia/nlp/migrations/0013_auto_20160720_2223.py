# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0012_auto_20160630_1231'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperengine',
            name='score_n_papers',
            field=models.PositiveIntegerField(default=10950),
        ),
    ]
