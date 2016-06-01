# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('nlp', '0016_threadengine_user_boost'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserFingerprint',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True, serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('added_after', models.DateField(blank=True, null=True)),
                ('embedding', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True)),
                ('authors_ids', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True)),
                ('authors_counts', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True)),
                ('journals_ids', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True)),
                ('journals_counts', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True)),
                ('user', models.OneToOneField(related_name='fingerprint', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
