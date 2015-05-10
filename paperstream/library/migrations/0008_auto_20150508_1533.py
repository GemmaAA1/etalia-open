# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0007_auto_20150507_2114'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='id_oth',
            field=models.CharField(null=True, db_index=True, max_length=32, blank=True, unique=True),
        ),
    ]
