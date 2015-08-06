# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0002_auto_20150708_0630'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mods',
            name='mods_path',
            field=models.CharField(max_length=256, default='/Users/nicolaspannetier/Google Drive/Projects/paperstream/paperstream_project/paperstream/nlp/mods'),
        ),
    ]
