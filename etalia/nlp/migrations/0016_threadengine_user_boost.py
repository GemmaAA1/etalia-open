# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0015_auto_20160531_0448'),
    ]

    operations = [
        migrations.AddField(
            model_name='threadengine',
            name='user_boost',
            field=models.FloatField(null=True, blank=True, default=0.05),
        ),
    ]
