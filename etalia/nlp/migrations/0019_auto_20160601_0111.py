# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0018_userfingerprint_model'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userfingerprint',
            name='authors_counts',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), null=True, size=300),
        ),
        migrations.AlterField(
            model_name='userfingerprint',
            name='authors_ids',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), null=True, size=300),
        ),
        migrations.AlterField(
            model_name='userfingerprint',
            name='journals_counts',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), null=True, size=300),
        ),
        migrations.AlterField(
            model_name='userfingerprint',
            name='journals_ids',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), null=True, size=300),
        ),
    ]
