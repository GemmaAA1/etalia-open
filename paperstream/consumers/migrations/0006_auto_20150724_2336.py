# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0005_auto_20150618_2215'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumer',
            name='day0',
            field=models.IntegerField(default=365),
        ),
    ]
