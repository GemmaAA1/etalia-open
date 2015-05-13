# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0016_auto_20150513_0501'),
    ]

    operations = [
        migrations.AddField(
            model_name='paper',
            name='id_pii',
            field=core.models.NullableCharField(verbose_name='Publisher ID', max_length=32, null=True, blank=True, unique=True, default=''),
        ),
    ]
