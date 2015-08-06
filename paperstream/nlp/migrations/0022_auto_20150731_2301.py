# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0021_auto_20150730_0021'),
    ]

    operations = [
        migrations.RenameField(
            model_name='model',
            old_name='alpha',
            new_name='_alpha',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='dbow_words',
            new_name='_dbow_words',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='dm',
            new_name='_dm',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='dm_concat',
            new_name='_dm_concat',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='dm_mean',
            new_name='_dm_mean',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='dm_tag_count',
            new_name='_dm_tag_count',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='hs',
            new_name='_hs',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='max_vocab_size',
            new_name='_max_vocab_size',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='min_count',
            new_name='_min_count',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='negative',
            new_name='_negative',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='sample',
            new_name='_sample',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='seed',
            new_name='_seed',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='size',
            new_name='_size',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='window',
            new_name='_window',
        ),
        migrations.RenameField(
            model_name='model',
            old_name='workers',
            new_name='_workers',
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='paper',
            field=models.ForeignKey(related_name='neighbors', to='library.Paper'),
        ),
    ]
