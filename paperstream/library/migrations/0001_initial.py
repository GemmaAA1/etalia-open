# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('first_name', models.CharField(default='', max_length=100, blank=True)),
                ('last_name', models.CharField(default='', max_length=100, blank=True)),
                ('email', models.EmailField(default='', max_length=254, blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='AuthorPosition',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('position', models.IntegerField(null=True, blank=True)),
                ('author', models.ForeignKey(to='library.Author')),
            ],
            options={
                'ordering': ['position'],
            },
        ),
        migrations.CreateModel(
            name='Journal',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('id_issn', models.CharField(default='', db_index=True, max_length=9, blank=True)),
                ('id_eissn', models.CharField(default='', db_index=True, max_length=9, blank=True)),
                ('id_arx', models.CharField(default='', db_index=True, max_length=32, blank=True)),
                ('id_oth', models.CharField(default='', db_index=True, max_length=32, blank=True)),
                ('title', models.CharField(default='', max_length=200, blank=True)),
                ('short_title', models.CharField(default='', max_length=100, blank=True)),
                ('url', models.URLField(default='', blank=True, null=True)),
                ('scope', models.TextField(default='', max_length=1000, blank=True)),
                ('language', models.CharField(default='', max_length=200, blank=True)),
                ('period', models.CharField(default='IRR', choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular')], max_length=200)),
                ('lib_size', models.IntegerField(default=0)),
                ('is_valid', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['title'],
            },
        ),
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('id_doi', models.CharField(default='', db_index=True, max_length=32, blank=True)),
                ('id_arx', models.CharField(default='', db_index=True, max_length=32, blank=True)),
                ('id_pmi', models.CharField(default='', db_index=True, max_length=32, blank=True)),
                ('id_oth', models.CharField(default='', db_index=True, max_length=32, blank=True)),
                ('title', models.CharField(default='', max_length=500, blank=True)),
                ('abstract', models.TextField(default='', blank=True)),
                ('volume', models.CharField(default='', max_length=200, blank=True)),
                ('issue', models.CharField(default='', max_length=200, blank=True)),
                ('page', models.CharField(default='', max_length=200, blank=True)),
                ('date', models.DateField(null=True, blank=True)),
                ('url', models.URLField(default='', blank=True, null=True)),
                ('is_aip', models.BooleanField(default=False)),
                ('is_pre_print', models.BooleanField(default=False)),
                ('is_valid', models.BooleanField(default=False)),
                ('authors', models.ManyToManyField(to='library.Author', through='library.AuthorPosition')),
                ('journal', models.ForeignKey(to='library.Journal', null=True, blank=True)),
            ],
            options={
                'ordering': ['-date'],
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(max_length=200, unique=True)),
                ('url', models.URLField()),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(to='library.Publisher', null=True, blank=True, default=''),
        ),
        migrations.AddField(
            model_name='authorposition',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AlterUniqueTogether(
            name='author',
            unique_together=set([('first_name', 'last_name')]),
        ),
    ]
