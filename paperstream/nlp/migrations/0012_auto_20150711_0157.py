# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150710_0053'),
        ('nlp', '0011_auto_20150709_2119'),
    ]

    operations = [
        migrations.CreateModel(
            name='JournalVectors',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=None, null=True)),
                ('journal', models.ForeignKey(to='library.Journal')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='model',
            name='status',
        ),
        migrations.AlterField(
            model_name='model',
            name='hs',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='journalvectors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
    ]
