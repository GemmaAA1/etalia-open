# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumer',
            name='day0',
            field=models.IntegerField(default=30),
        ),
    ]
