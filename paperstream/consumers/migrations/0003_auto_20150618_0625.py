# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0002_auto_20150521_2314'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumer',
            name='day0',
            field=models.IntegerField(default=180),
        ),
    ]
