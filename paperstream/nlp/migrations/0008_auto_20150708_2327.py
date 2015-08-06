# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0007_auto_20150708_2300'),
    ]

    operations = [
        migrations.AlterField(
            model_name='papervectors',
            name='vector',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=None, null=True),
        ),
    ]
