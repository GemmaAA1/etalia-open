# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models
import library.validators


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_auto_20150510_0859'),
    ]

    operations = [
        migrations.AlterField(
            model_name='author',
            name='last_name',
            field=models.CharField(max_length=100, validators=[library.validators.validate_author_names]),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_eissn',
            field=core.models.NullableCharField(default=None, validators=[library.validators.validate_issn], max_length=9, blank=True, null=True, unique=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_issn',
            field=core.models.NullableCharField(default=None, validators=[library.validators.validate_issn], max_length=9, blank=True, null=True, unique=True),
        ),
    ]
