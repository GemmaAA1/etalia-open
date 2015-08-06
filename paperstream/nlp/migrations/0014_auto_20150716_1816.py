# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0013_auto_20150716_0517'),
    ]

    operations = [
        migrations.AlterField(
            model_name='model',
            name='size',
            field=models.IntegerField(default=100, validators=[django.core.validators.MaxValueValidator(300)]),
        ),
        migrations.AlterField(
            model_name='papervectors',
            name='paper',
            field=models.ForeignKey(to='library.Paper', related_name='vectors'),
        ),
        migrations.AlterField(
            model_name='papervectors',
            name='vector',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300),
        ),
    ]
