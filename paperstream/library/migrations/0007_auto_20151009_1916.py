# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20151009_0629'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='is_trusted',
            field=models.BooleanField(default=False, db_index=True),
        ),
    ]
