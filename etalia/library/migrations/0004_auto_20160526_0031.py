# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0010_auto_20160526_0031'),
        ('library', '0003_journal_lib_size'),
    ]

    operations = [
        migrations.CreateModel(
            name='JournalNeighbors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), blank=True, null=True, size=10)),
                ('journal', models.ForeignKey(to='library.Journal', related_name='neighbors')),
                ('ms', models.ForeignKey(to='nlp.MostSimilar')),
            ],
        ),
        migrations.CreateModel(
            name='JournalVectors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300, db_index=True)),
                ('journal', models.ForeignKey(to='library.Journal', related_name='vectors')),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
        ),
        migrations.CreateModel(
            name='PaperNeighbors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(verbose_name='Days from right now', default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')])),
                ('neighbors', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), blank=True, null=True, size=10)),
                ('ms', models.ForeignKey(to='nlp.MostSimilar')),
                ('paper', models.ForeignKey(to='library.Paper', related_name='neighbors')),
            ],
        ),
        migrations.CreateModel(
            name='PaperVectors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300, db_index=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('paper', models.ForeignKey(to='library.Paper', related_name='vectors')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([('paper', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('time_lapse', 'paper', 'ms')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([('journal', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([('journal', 'ms')]),
        ),
    ]
