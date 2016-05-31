# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields
import etalia.nlp.mixins


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0029_auto_20160513_0042'),
        ('library', '0003_journal_lib_size'),
        ('nlp', '0010_auto_20160526_0031'),
    ]

    operations = [
        migrations.CreateModel(
            name='JournalNeighbors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, size=10, null=True, base_field=models.IntegerField(null=True))),
                ('journal', models.ForeignKey(related_name='neighbors', to='library.Journal')),
            ],
        ),
        migrations.CreateModel(
            name='JournalVectors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(size=300, null=True, base_field=models.FloatField(null=True), db_index=True)),
                ('journal', models.ForeignKey(related_name='vectors', to='library.Journal')),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
        ),
        migrations.CreateModel(
            name='PaperEngine',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('upload_state', models.CharField(default='IDL ', max_length=3, choices=[('IDL', 'Idle'), ('ING', 'Uploading')])),
                ('journal_ratio', models.FloatField(default=0.0, choices=[('0.0', '0 %'), ('0.05', '5 %'), ('0.10', '10 %'), ('0.15', '15 %'), ('0.20', '20 %'), ('0.25', '25 %'), ('0.30', '30 %')])),
                ('is_active', models.BooleanField(default=False)),
                ('model', models.ForeignKey(related_name='pe', to='nlp.Model')),
            ],
            bases=(models.Model, etalia.nlp.mixins.S3Mixin),
        ),
        migrations.CreateModel(
            name='PaperNeighbors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Days from right now')),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, size=10, null=True, base_field=models.IntegerField(null=True))),
                ('pe', models.ForeignKey(to='nlp.PaperEngine')),
                ('paper', models.ForeignKey(related_name='neighbors', to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='PaperVectors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(size=300, null=True, base_field=models.FloatField(null=True), db_index=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('paper', models.ForeignKey(related_name='vectors', to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='ThreadEngine',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('upload_state', models.CharField(default='IDL ', max_length=3, choices=[('IDL', 'Idle'), ('ING', 'Uploading')])),
                ('is_active', models.BooleanField(default=False)),
                ('model', models.ForeignKey(related_name='ms_thread', to='nlp.Model')),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, etalia.nlp.mixins.S3Mixin),
        ),
        migrations.CreateModel(
            name='ThreadNeighbors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Days from right now')),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, size=10, null=True, base_field=models.IntegerField(null=True))),
                ('pe', models.ForeignKey(to='nlp.ThreadEngine')),
                ('thread', models.ForeignKey(to='threads.Thread')),
            ],
        ),
        migrations.CreateModel(
            name='ThreadVectors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(size=300, null=True, base_field=models.FloatField(null=True), db_index=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('thread', models.ForeignKey(related_name='vectors', to='threads.Thread')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='mostsimilar',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='mostsimilar',
            name='model',
        ),
        migrations.RemoveField(
            model_name='mostsimilarthread',
            name='model',
        ),
        migrations.DeleteModel(
            name='MostSimilar',
        ),
        migrations.DeleteModel(
            name='MostSimilarThread',
        ),
        migrations.AddField(
            model_name='journalneighbors',
            name='pe',
            field=models.ForeignKey(to='nlp.PaperEngine'),
        ),
        migrations.AlterUniqueTogether(
            name='threadvectors',
            unique_together=set([('thread', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadneighbors',
            unique_together=set([('time_lapse', 'thread', 'pe')]),
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([('paper', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('time_lapse', 'paper', 'pe')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperengine',
            unique_together=set([('model', 'journal_ratio')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([('journal', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([('journal', 'pe')]),
        ),
    ]
