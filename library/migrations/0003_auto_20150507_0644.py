# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_auto_20150507_0033'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='identifiers',
            field=jsonfield.fields.JSONField(unique=True),
        ),
    ]
