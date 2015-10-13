# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0009_auto_20151013_1958'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='mostsimilarstatus',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='mostsimilarstatus',
            name='ms',
        ),
        migrations.DeleteModel(
            name='MostSimilarStatus',
        ),
    ]
