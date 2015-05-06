# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Creator',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('first_name', models.CharField(blank=True, default='', max_length=100)),
                ('last_name', models.CharField(blank=True, default='', max_length=100)),
                ('email', models.EmailField(blank=True, null=True, max_length=254)),
            ],
        ),
        migrations.CreateModel(
            name='CreatorPosition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('position', models.IntegerField(blank=True, null=True)),
                ('creator', models.ForeignKey(to='library.Creator')),
            ],
            options={
                'ordering': ['position'],
            },
        ),
        migrations.CreateModel(
            name='Journal',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(max_length=200)),
                ('short_title', models.CharField(blank=True, null=True, default='', max_length=100)),
                ('issn', models.CharField(blank=True, null=True, default='', max_length=10)),
                ('e_issn', models.CharField(blank=True, null=True, default='', max_length=10)),
                ('ext_id', models.CharField(blank=True, null=True, default='', max_length=30)),
                ('url', models.URLField(blank=True, null=True)),
                ('scope', models.TextField(blank=True, default='', max_length=1000)),
                ('language', models.CharField(blank=True, null=True, default='', max_length=200)),
                ('period', models.CharField(choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular')], default='IRR', max_length=200)),
                ('paper_counter', models.IntegerField(default=0)),
                ('is_valid', models.BooleanField(default=False)),
            ],
        ),
        migrations.CreateModel(
            name='NewPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('lib_size', models.IntegerField(default=0)),
            ],
        ),
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(blank=True, max_length=500)),
                ('abstract', models.TextField(blank=True)),
                ('volume', models.CharField(blank=True, max_length=200)),
                ('issue', models.CharField(blank=True, max_length=200)),
                ('page', models.CharField(blank=True, max_length=200)),
                ('date', models.DateField(blank=True, null=True)),
                ('doi', models.CharField(blank=True, max_length=200)),
                ('ext_id', models.CharField(blank=True, max_length=200)),
                ('url', models.URLField(blank=True, null=True)),
                ('date_added', models.DateTimeField(auto_now_add=True)),
                ('is_aip', models.BooleanField(default=False)),
                ('is_pre_print', models.BooleanField(default=False)),
                ('is_valid', models.BooleanField(default=False)),
                ('creators', models.ManyToManyField(through='library.CreatorPosition', to='library.Creator')),
                ('journal', models.ForeignKey(null=True, to='library.Journal')),
            ],
            options={
                'ordering': ['-date'],
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
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
            field=models.ForeignKey(null=True, to='library.Publisher'),
        ),
        migrations.AddField(
            model_name='creatorposition',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AlterUniqueTogether(
            name='creator',
            unique_together=set([('first_name', 'last_name')]),
        ),
        migrations.AlterUniqueTogether(
            name='paper',
            unique_together=set([('doi', 'ext_id')]),
        ),
        migrations.AlterUniqueTogether(
            name='journal',
            unique_together=set([('issn', 'ext_id')]),
        ),
    ]
