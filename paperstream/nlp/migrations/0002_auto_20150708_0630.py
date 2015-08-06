# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mods',
            name='max_vocab_size',
            field=models.IntegerField(default=None, null=True, blank=True),
        ),
    ]
