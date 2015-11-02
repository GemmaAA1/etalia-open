# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Stream',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=100)),
                ('state', models.CharField(default='NON', blank=True, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3)),
                ('message', models.CharField(default='', blank=True, max_length=127)),
            ],
        ),
        migrations.CreateModel(
            name='StreamMatches',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('is_score_computed', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.CreateModel(
            name='StreamSeeds',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
            ],
        ),
        migrations.CreateModel(
            name='StreamVector',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300)),
            ],
        ),
        migrations.CreateModel(
            name='Trend',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=100)),
                ('top_n_closest', models.IntegerField(default=5000)),
                ('time_lapse_top_altmetric', models.IntegerField(default=15)),
                ('score_threshold', models.FloatField(default=1000.0)),
            ],
        ),
        migrations.CreateModel(
            name='TrendMatches',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('trend_feed', models.ForeignKey(to='feeds.Trend')),
            ],
        ),
        migrations.AddField(
            model_name='trend',
            name='matches',
            field=models.ManyToManyField(through='feeds.TrendMatches', to='library.Paper'),
        ),
    ]
