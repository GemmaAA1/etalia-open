# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0009_auto_20150709_1537'),
    ]

    operations = [
        migrations.AlterField(
            model_name='model',
            name='doc2vec_path',
            field=models.CharField(max_length=256, default='/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/nlp/mods'),
        ),
    ]
