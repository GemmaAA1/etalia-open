# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import library.validators


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_auto_20150510_0114'),
    ]

    operations = [
        migrations.AlterField(
            model_name='author',
            name='last_name',
            field=models.CharField(max_length=100),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_arx',
            field=models.CharField(default='', blank=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_eissn',
            field=models.CharField(default='', blank=True, max_length=9, validators=[library.validators.validate_id_eissn]),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_issn',
            field=models.CharField(default='', blank=True, max_length=9, validators=[library.validators.validate_id_issn]),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_oth',
            field=models.CharField(default='', blank=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='journal',
            name='short_title',
            field=models.CharField(blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='journal',
            name='title',
            field=models.CharField(max_length=200),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(default='', blank=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(default='', blank=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(default='', blank=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(default='', blank=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='title',
            field=models.CharField(max_length=500),
        ),
    ]
