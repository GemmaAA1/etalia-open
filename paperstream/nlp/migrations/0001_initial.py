# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_paper_id_isbn'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fingerprint',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vec', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(), size=None)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Mods',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(unique=True, max_length=128)),
                ('mods_path', models.CharField(max_length=256, default='/Users/nicolaspannetier/Google Drive/Projects/paperstream/paperstream_project/paperstream/nlp/data')),
                ('status', models.CharField(choices=[('UNT', 'Untrained'), ('VOC', 'Building Vocabulary'), ('TRA', 'Training'), ('SAV', 'Saving'), ('LOA', 'Loading'), ('IDL', 'Idle')], max_length=3, default='UNT')),
                ('dm', models.IntegerField(default=1)),
                ('size', models.IntegerField(default=100)),
                ('window', models.IntegerField(default=8)),
                ('alpha', models.FloatField(default=0.025)),
                ('seed', models.IntegerField(default=0)),
                ('min_count', models.IntegerField(default=5)),
                ('max_vocab_size', models.IntegerField(default=None)),
                ('sample', models.FloatField(default=0.0)),
                ('workers', models.IntegerField(default=1)),
                ('hs', models.IntegerField(default=1)),
                ('negative', models.IntegerField(default=0)),
                ('dm_mean', models.IntegerField(default=0)),
                ('dm_concat', models.IntegerField(default=0)),
                ('dm_tag_count', models.IntegerField(default=1)),
                ('dbow_words', models.IntegerField(default=0)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='fingerprint',
            name='mod',
            field=models.ForeignKey(to='nlp.Mods'),
        ),
        migrations.AddField(
            model_name='fingerprint',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
    ]
