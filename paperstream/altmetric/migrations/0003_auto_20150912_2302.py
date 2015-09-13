# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('altmetric', '0002_auto_20150909_1624'),
    ]

    operations = [
        migrations.AlterField(
            model_name='altmetricmodel',
            name='altmetric_jid',
            field=models.CharField(max_length=32, blank=True),
        ),
        migrations.AlterField(
            model_name='altmetricmodel',
            name='score',
            field=models.FloatField(default=0.0, db_index=True),
        ),
    ]
