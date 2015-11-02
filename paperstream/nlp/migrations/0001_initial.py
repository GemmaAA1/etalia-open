# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators
import paperstream.nlp.mixins
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='JournalNeighbors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, base_field=models.IntegerField(null=True), null=True, size=10)),
                ('journal', models.ForeignKey(related_name='neighbors', to='library.Journal')),
            ],
        ),
        migrations.CreateModel(
            name='JournalVectors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300, db_index=True)),
                ('journal', models.ForeignKey(related_name='vectors', to='library.Journal')),
            ],
        ),
        migrations.CreateModel(
            name='Model',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(max_length=128, unique=True)),
                ('is_active', models.BooleanField(default=False)),
                ('upload_state', models.CharField(default='IDL', choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3)),
                ('_dm', models.IntegerField(default=1)),
                ('_size', models.IntegerField(default=128, validators=[django.core.validators.MaxValueValidator(300)])),
                ('_window', models.IntegerField(default=8)),
                ('_alpha', models.FloatField(default=0.025)),
                ('_seed', models.IntegerField(default=0)),
                ('_min_count', models.IntegerField(default=2)),
                ('_max_vocab_size', models.IntegerField(default=None, blank=True, null=True)),
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
            bases=(models.Model, paperstream.nlp.mixins.S3Mixin),
        ),
        migrations.CreateModel(
            name='MostSimilar',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('upload_state', models.CharField(default='IDL ', choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3)),
                ('journal_ratio', models.FloatField(default=0.0)),
                ('is_active', models.BooleanField(default=False)),
                ('model', models.ForeignKey(related_name='ms', to='nlp.Model')),
            ],
            bases=(models.Model, paperstream.nlp.mixins.S3Mixin),
        ),
        migrations.CreateModel(
            name='PaperNeighbors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (365, 'Year'), (-1, 'All')], verbose_name='Days from right now')),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, base_field=models.IntegerField(null=True), null=True, size=10)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('paper', models.ForeignKey(related_name='neighbors', to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='PaperVectors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300, db_index=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('paper', models.ForeignKey(related_name='vectors', to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='TextField',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('text_field', models.CharField(choices=[('title', 'Title'), ('abstract', 'Abstract')], max_length=100)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='model',
            name='text_fields',
            field=models.ManyToManyField(related_name='text_fields', to='nlp.TextField'),
        ),
        migrations.AddField(
            model_name='journalvectors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='journalneighbors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([('paper', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('time_lapse', 'paper', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='mostsimilar',
            unique_together=set([('model', 'journal_ratio')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([('journal', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([('journal', 'model')]),
        ),
    ]
