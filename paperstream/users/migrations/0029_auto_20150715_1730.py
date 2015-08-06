# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0028_auto_20150710_0045'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userlibpaper',
            name='date_created',
            field=models.DateField(default=None),
        ),
        migrations.AlterField(
            model_name='userlibpaper',
            name='date_last_modified',
            field=models.DateField(default=None),
        ),
    ]
