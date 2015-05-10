# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import library.models
import library.validators


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20150510_0745'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='id_arx',
            field=library.models.NullableCharField(null=True, blank=True, unique=True, default=None, max_length=32),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_eissn',
            field=library.models.NullableCharField(null=True, validators=[library.validators.validate_id_eissn], unique=True, default=None, max_length=9, blank=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_issn',
            field=library.models.NullableCharField(null=True, validators=[library.validators.validate_id_issn], unique=True, default=None, max_length=9, blank=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_oth',
            field=library.models.NullableCharField(null=True, blank=True, unique=True, default=None, max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=library.models.NullableCharField(null=True, blank=True, unique=True, default='', max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=library.models.NullableCharField(null=True, blank=True, unique=True, default='', max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=library.models.NullableCharField(null=True, blank=True, unique=True, default='', max_length=32),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=library.models.NullableCharField(null=True, blank=True, unique=True, default='', max_length=32),
        ),
    ]
