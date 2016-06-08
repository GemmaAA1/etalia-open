# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0003_paperengine_score_altmetric_boost'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfingerprint',
            name='thread_embedding',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300),
        ),
        migrations.AddField(
            model_name='userfingerprint',
            name='threads_users_counts',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), null=True, size=300),
        ),
        migrations.AddField(
            model_name='userfingerprint',
            name='threads_users_ids',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), null=True, size=300),
        ),
    ]
