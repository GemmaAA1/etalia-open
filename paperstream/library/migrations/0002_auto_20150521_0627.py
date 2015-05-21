# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=core.models.NullableCharField(verbose_name='Arxiv', null=True, max_length=64, default='', unique=True, blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=core.models.NullableCharField(verbose_name='DOI', null=True, max_length=64, default='', unique=True, blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=core.models.NullableCharField(verbose_name='Other ID', null=True, max_length=64, default='', unique=True, blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pii',
            field=core.models.NullableCharField(verbose_name='PII', null=True, default='', blank=True, max_length=64),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=core.models.NullableCharField(verbose_name='PMID', null=True, max_length=64, default='', unique=True, blank=True),
        ),
    ]
