# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('first_name', models.CharField(max_length=100, blank=True, default='')),
                ('last_name', models.CharField(max_length=100, blank=True, default='')),
                ('email', models.EmailField(null=True, max_length=254, blank=True, default='')),
            ],
        ),
        migrations.CreateModel(
            name='AuthorPosition',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
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
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('id_key', models.CharField(max_length=200, choices=[('ISSN', 'ISSN'), ('ARX', 'Arxiv'), ('OTH', 'Other')], default='OTH', db_index=True)),
                ('id_val', models.CharField(max_length=240, db_index=True)),
                ('title', models.CharField(max_length=200, blank=True, default='')),
                ('short_title', models.CharField(max_length=100, blank=True, default='')),
                ('issn', models.CharField(max_length=10, blank=True, default='')),
                ('e_issn', models.CharField(max_length=10, blank=True, default='')),
                ('ext_id', models.CharField(max_length=30, blank=True, default='')),
                ('url', models.URLField(null=True, blank=True, default='')),
                ('scope', models.TextField(max_length=1000, blank=True, default='')),
                ('language', models.CharField(max_length=200, blank=True, default='')),
                ('period', models.CharField(max_length=200, choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular')], default='IRR')),
                ('lib_size', models.IntegerField(default=0)),
                ('is_valid', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['title'],
            },
        ),
        migrations.CreateModel(
            name='NewPaper',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('lib_size', models.IntegerField(default=0)),
            ],
        ),
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('identifiers', jsonfield.fields.JSONField(default='o')),
                ('title', models.CharField(max_length=500, blank=True, default='')),
                ('abstract', models.TextField(blank=True, default='')),
                ('volume', models.CharField(max_length=200, blank=True, default='')),
                ('issue', models.CharField(max_length=200, blank=True, default='')),
                ('page', models.CharField(max_length=200, blank=True, default='')),
                ('date', models.DateField(null=True, blank=True)),
                ('url', models.URLField(null=True, blank=True, default='')),
                ('date_added', models.DateTimeField(auto_now_add=True)),
                ('is_aip', models.BooleanField(default=False)),
                ('is_pre_print', models.BooleanField(default=False)),
                ('is_valid', models.BooleanField(default=False)),
                ('authors', models.ManyToManyField(through='library.AuthorPosition', to='library.Author')),
                ('journal', models.ForeignKey(to='library.Journal', null=True, blank=True)),
            ],
            options={
                'ordering': ['-date'],
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('name', models.CharField(unique=True, max_length=200)),
                ('url', models.URLField()),
            ],
        ),
        migrations.AddField(
            model_name='newpaper',
            name='paper',
            field=models.OneToOneField(to='library.Paper'),
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
        migrations.AlterUniqueTogether(
            name='journal',
            unique_together=set([('id_key', 'id_val')]),
        ),
    ]
