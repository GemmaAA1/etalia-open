# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0005_auto_20150710_0045'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='title',
            field=models.CharField(db_index=True, default='', max_length=500),
        ),
    ]
