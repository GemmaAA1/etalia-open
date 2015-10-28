# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0012_auto_20151027_0521'),
    ]

    operations = [
        migrations.AddField(
            model_name='mostsimilar',
            name='journal_ratio',
            field=models.FloatField(default=0.0),
        ),
    ]
