# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='identifiers',
            field=jsonfield.fields.JSONField(db_index=True, unique=True),
        ),
    ]
