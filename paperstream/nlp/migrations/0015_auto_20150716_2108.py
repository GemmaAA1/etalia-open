# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0014_auto_20150716_1816'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journalvectors',
            name='vector',
            field=django.contrib.postgres.fields.ArrayField(null=True, base_field=models.FloatField(null=True), size=300),
        ),
    ]
