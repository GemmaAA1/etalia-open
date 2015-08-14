# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0031_auto_20150806_0308'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journalvectors',
            name='journal',
            field=models.ForeignKey(related_name='vectors', to='library.Journal'),
        ),
    ]
