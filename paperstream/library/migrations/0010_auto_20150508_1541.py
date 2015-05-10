# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0009_auto_20150508_1537'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='journal',
            name='e_issn',
        ),
        migrations.RemoveField(
            model_name='journal',
            name='ext_id',
        ),
        migrations.RemoveField(
            model_name='journal',
            name='issn',
        ),
        migrations.AddField(
            model_name='journal',
            name='id_eissn',
            field=models.CharField(blank=True, null=True, max_length=9, db_index=True, unique=True),
        ),
    ]
