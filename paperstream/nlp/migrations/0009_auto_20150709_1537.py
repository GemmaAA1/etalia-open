# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0008_auto_20150708_2327'),
    ]

    operations = [
        migrations.AddField(
            model_name='model',
            name='is_active',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='model',
            name='status',
            field=models.CharField(choices=[('UNT', 'Untrained'), ('VOC', 'Building Vocabulary'), ('TRA', 'Training'), ('SAV', 'Saving'), ('LOA', 'Loading'), ('POP', 'Populating'), ('IDL', 'Idle'), ('USE', 'Usable')], max_length=3, default='UNT'),
        ),
    ]
