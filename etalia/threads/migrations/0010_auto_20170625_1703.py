# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0009_auto_20160922_2029'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pubpeer',
            name='doi',
            field=models.CharField(blank=True, max_length=64, verbose_name='DOI', default='', db_index=True),
        ),
    ]
