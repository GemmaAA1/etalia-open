# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0004_auto_20150708_0714'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='Mods',
            new_name='Model',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='mods_path',
            new_name='doc2vec_path',
        ),
    ]
