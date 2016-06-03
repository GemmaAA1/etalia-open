# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Stream',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=128)),
                ('state', models.CharField(blank=True, default='NON', choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3)),
                ('last_update', models.DateTimeField(blank=True, default=None, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='StreamPapers',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('date', models.DateField()),
                ('new', models.BooleanField(default=True)),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.CreateModel(
            name='Trend',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=128)),
                ('state', models.CharField(blank=True, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3)),
                ('last_update', models.DateTimeField(blank=True, default=None, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='TrendPapers',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('date', models.DateField()),
                ('new', models.BooleanField(default=True)),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('trend', models.ForeignKey(to='feeds.Trend')),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.AddField(
            model_name='trend',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='feeds.TrendPapers'),
        ),
    ]
