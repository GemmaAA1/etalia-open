# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='JournalVectors',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(null=True, base_field=models.FloatField(null=True), db_index=True, size=300)),
                ('journal', models.ForeignKey(to='library.Journal', related_name='vectors')),
            ],
        ),
        migrations.CreateModel(
            name='LSH',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('state', models.CharField(choices=[('NON', 'None'), ('BUS', 'Busy'), ('IDL', 'Idle')], max_length=3, default='NON')),
                ('time_lapse', models.IntegerField(choices=[(61, '2 Months'), (-1, 'All')], default=-1, verbose_name='Days from right now')),
            ],
        ),
        migrations.CreateModel(
            name='Model',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(unique=True, max_length=128)),
                ('is_active', models.BooleanField(default=False)),
                ('_dm', models.IntegerField(default=1)),
                ('_size', models.IntegerField(default=128, validators=[django.core.validators.MaxValueValidator(300)])),
                ('_window', models.IntegerField(default=8)),
                ('_alpha', models.FloatField(default=0.025)),
                ('_seed', models.IntegerField(default=0)),
                ('_min_count', models.IntegerField(default=2)),
                ('_max_vocab_size', models.IntegerField(null=True, blank=True, default=None)),
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
        ),
        migrations.CreateModel(
            name='PaperNeighbors',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(null=True, blank=True, base_field=models.IntegerField(null=True), size=10)),
                ('lsh', models.ForeignKey(to='nlp.LSH')),
                ('paper', models.ForeignKey(to='library.Paper', related_name='neighbors')),
            ],
        ),
        migrations.CreateModel(
            name='PaperVectors',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(null=True, base_field=models.FloatField(null=True), db_index=True, size=300)),
                ('is_in_full_lsh', models.BooleanField(default=False)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('paper', models.ForeignKey(to='library.Paper', related_name='vectors')),
            ],
        ),
        migrations.CreateModel(
            name='TextField',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
            field=models.ManyToManyField(to='nlp.TextField', related_name='text_fields'),
        ),
        migrations.AddField(
            model_name='lsh',
            name='model',
            field=models.ForeignKey(to='nlp.Model', related_name='lsh'),
        ),
        migrations.AddField(
            model_name='journalvectors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([('paper', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('lsh', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='lsh',
            unique_together=set([('model', 'time_lapse')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([('journal', 'model')]),
        ),
    ]
