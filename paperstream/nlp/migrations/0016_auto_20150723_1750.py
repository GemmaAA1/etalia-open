# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0015_auto_20150716_2108'),
    ]

    operations = [
        migrations.CreateModel(
            name='LSHFModel',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('state', models.CharField(default='', max_length=3, choices=[('NON', 'None'), ('BUS', 'Busy'), ('IDL', 'Idle')])),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='model',
            name='data_path',
        ),
        migrations.RemoveField(
            model_name='model',
            name='doc2vec_path',
        ),
        migrations.AddField(
            model_name='papervectors',
            name='is_in_lshf',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='journalvectors',
            name='vector',
            field=django.contrib.postgres.fields.ArrayField(size=300, base_field=models.FloatField(null=True), db_index=True, null=True),
        ),
        migrations.AlterField(
            model_name='papervectors',
            name='vector',
            field=django.contrib.postgres.fields.ArrayField(size=300, base_field=models.FloatField(null=True), db_index=True, null=True),
        ),
        migrations.AddField(
            model_name='lshfmodel',
            name='model',
            field=models.ForeignKey(to='nlp.Model', unique=True),
        ),
    ]
