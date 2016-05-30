# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0013_auto_20160526_2359'),
    ]

    operations = [
        migrations.AddField(
            model_name='paperengine',
            name='journal_boost',
            field=models.FloatField(default=0.05, blank=True, null=True),
        ),
    ]
