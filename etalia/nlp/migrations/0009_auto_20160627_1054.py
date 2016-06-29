# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0008_auto_20160627_1049'),
    ]

    operations = [
        migrations.AlterField(
            model_name='threadengine',
            name='score_paper_boost',
            field=models.FloatField(blank=True, null=True, default=0.4),
        ),
        migrations.AlterField(
            model_name='threadengine',
            name='score_user_boost',
            field=models.FloatField(blank=True, null=True, default=0.2),
        ),
    ]
