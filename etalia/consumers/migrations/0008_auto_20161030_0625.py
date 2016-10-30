# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0007_auto_20161013_1829'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumer',
            name='day0',
            field=models.IntegerField(default=7),
        ),
    ]
