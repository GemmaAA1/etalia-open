# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20150528_0034'),
    ]

    operations = [
        migrations.AddField(
            model_name='paper',
            name='id_isbn',
            field=core.models.NullableCharField(db_index=True, unique=True, max_length=64, verbose_name='ISBN', blank=True, null=True, default=''),
        ),
    ]
