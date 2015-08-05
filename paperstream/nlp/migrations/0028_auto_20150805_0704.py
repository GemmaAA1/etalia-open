# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0027_auto_20150804_1717'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lsh',
            name='state',
            field=models.CharField(max_length=3, choices=[('NON', 'None'), ('BUS', 'Busy'), ('IDL', 'Idle')], default='NON'),
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='neighbors',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), null=True, size=10, blank=True),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('lsh', 'paper')]),
        ),
    ]
