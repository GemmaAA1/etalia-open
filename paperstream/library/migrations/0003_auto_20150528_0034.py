# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models
import library.validators


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_auto_20150527_2349'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='id_arx',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', blank=True, verbose_name='Arxiv ID', null=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_eissn',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', validators=[library.validators.validate_issn], blank=True, verbose_name='e-ISSN', null=True, max_length=9),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_issn',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', validators=[library.validators.validate_issn], blank=True, verbose_name='ISSN', null=True, max_length=9),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_oth',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', blank=True, verbose_name='Other ID', null=True, max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', blank=True, verbose_name='Arxiv', null=True, max_length=64),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', blank=True, verbose_name='DOI', null=True, max_length=64),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', blank=True, verbose_name='Other ID', null=True, max_length=64),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pii',
            field=core.models.NullableCharField(db_index=True, default='', blank=True, verbose_name='PII', null=True, max_length=64),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=core.models.NullableCharField(unique=True, db_index=True, default='', blank=True, verbose_name='PMID', null=True, max_length=64),
        ),
    ]
