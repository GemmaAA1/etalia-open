# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import etalia.nlp.mixins
import etalia.nlp.models.library
import django.contrib.postgres.fields
import django.core.validators
import etalia.nlp.models.threads


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='JournalNeighbors',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, base_field=models.IntegerField(null=True), size=10, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='JournalVectors',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), db_index=True, size=300, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Model',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(unique=True, max_length=128)),
                ('is_active', models.BooleanField(default=False)),
                ('upload_state', models.CharField(default='IDL', choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3)),
                ('_dm', models.IntegerField(default=1)),
                ('_size', models.IntegerField(default=128, validators=[django.core.validators.MaxValueValidator(300)])),
                ('_window', models.IntegerField(default=8)),
                ('_alpha', models.FloatField(default=0.025)),
                ('_seed', models.IntegerField(default=0)),
                ('_min_count', models.IntegerField(default=2)),
                ('_max_vocab_size', models.IntegerField(blank=True, default=None, null=True)),
                ('_sample', models.FloatField(default=0.0)),
                ('_workers', models.IntegerField(default=1)),
                ('_hs', models.IntegerField(default=1)),
                ('_negative', models.IntegerField(default=0)),
                ('_dm_mean', models.IntegerField(default=0)),
                ('_dm_concat', models.IntegerField(default=0)),
                ('_dm_tag_count', models.IntegerField(default=1)),
                ('_dbow_words', models.IntegerField(default=0)),
            ],
            options={
                'ordering': ['name'],
            },
            bases=(etalia.nlp.models.threads.ModelThreadMixin, etalia.nlp.models.library.ModelLibraryMixin, etalia.nlp.mixins.S3Mixin, models.Model),
        ),
        migrations.CreateModel(
            name='PaperEngine',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('upload_state', models.CharField(default='IDL ', choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3)),
                ('is_active', models.BooleanField(default=False)),
                ('journal_boost', models.FloatField(blank=True, default=0.05, null=True)),
                ('score_author_boost', models.FloatField(default=0.05)),
                ('score_journal_boost', models.FloatField(default=0.05)),
            ],
            options={
                'abstract': False,
            },
            bases=(etalia.nlp.models.library.PaperEngineScoringMixin, etalia.nlp.mixins.S3Mixin, models.Model),
        ),
        migrations.CreateModel(
            name='PaperNeighbors',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Days from right now')),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, base_field=models.IntegerField(null=True), size=10, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='PaperVectors',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), db_index=True, size=300, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='TextField',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('text_field', models.CharField(choices=[('title', 'Title'), ('abstract', 'Abstract')], max_length=100)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadEngine',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('upload_state', models.CharField(default='IDL ', choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3)),
                ('is_active', models.BooleanField(default=False)),
                ('user_boost', models.FloatField(blank=True, default=0.05, null=True)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, etalia.nlp.mixins.S3Mixin),
        ),
        migrations.CreateModel(
            name='ThreadNeighbors',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Days from right now')),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, base_field=models.IntegerField(null=True), size=10, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='ThreadVectors',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), db_index=True, size=300, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserFingerprint',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=128)),
                ('added_after', models.DateField(blank=True, null=True)),
                ('embedding', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True)),
                ('authors_ids', django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), size=300, null=True)),
                ('authors_counts', django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), size=300, null=True)),
                ('journals_ids', django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), size=300, null=True)),
                ('journals_counts', django.contrib.postgres.fields.ArrayField(base_field=models.PositiveIntegerField(null=True), size=300, null=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
        ),
    ]
