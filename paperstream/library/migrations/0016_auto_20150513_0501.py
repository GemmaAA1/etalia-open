# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models
import library.validators


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0015_auto_20150511_1842'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='id_arx',
            field=core.models.NullableCharField(null=True, blank=True, default=None, max_length=32, unique=True, verbose_name='Arxiv ID'),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_eissn',
            field=core.models.NullableCharField(null=True, blank=True, default=None, max_length=9, validators=[library.validators.validate_issn], unique=True, verbose_name='e-ISSN'),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_issn',
            field=core.models.NullableCharField(null=True, blank=True, default=None, max_length=9, validators=[library.validators.validate_issn], unique=True, verbose_name='ISSN'),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_oth',
            field=core.models.NullableCharField(null=True, blank=True, default=None, max_length=32, unique=True, verbose_name='Other ID'),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=core.models.NullableCharField(null=True, blank=True, default='', max_length=32, unique=True, verbose_name='Arxiv'),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=core.models.NullableCharField(null=True, blank=True, default='', max_length=32, unique=True, verbose_name='DOI'),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=core.models.NullableCharField(null=True, blank=True, default='', max_length=32, unique=True, verbose_name='Other ID'),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=core.models.NullableCharField(null=True, blank=True, default='', max_length=32, unique=True, verbose_name='PMID'),
        ),
    ]
