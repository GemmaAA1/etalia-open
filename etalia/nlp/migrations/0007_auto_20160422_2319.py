# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import etalia.nlp.mixins
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0027_auto_20160421_1748'),
        ('nlp', '0006_auto_20160421_0440'),
    ]

    operations = [
        migrations.CreateModel(
            name='MostSimilarThread',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('upload_state', models.CharField(default='IDL ', choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3)),
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
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(verbose_name='Days from right now', choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], default=-1)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), size=10, null=True, blank=True)),
                ('thread', models.ForeignKey(to='threads.Thread')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='threadneighbor',
            name='thread',
        ),
        migrations.DeleteModel(
            name='ThreadNeighbor',
        ),
    ]
