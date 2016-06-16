# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0004_auto_20160608_1800'),
    ]

    operations = [
        migrations.AddField(
            model_name='threadengine',
            name='paper_boost',
            field=models.FloatField(null=True, default=0.2, blank=True),
        ),
        migrations.AddField(
            model_name='threadengine',
            name='score_thread_embedding_paper_weight',
            field=models.FloatField(null=True, default=1.0, blank=True),
        ),
        migrations.AlterField(
            model_name='threadengine',
            name='user_boost',
            field=models.FloatField(null=True, default=0.1, blank=True),
        ),
    ]
