# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0009_auto_20150511_1712'),
    ]

    operations = [
        migrations.AlterField(
            model_name='publisher',
            name='url',
            field=models.URLField(default=True, blank=True),
        ),
    ]
