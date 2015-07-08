# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0003_auto_20150708_0646'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mods',
            name='min_count',
            field=models.IntegerField(default=2),
        ),
    ]
