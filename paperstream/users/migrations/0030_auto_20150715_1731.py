# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0029_auto_20150715_1730'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userlibpaper',
            name='date_created',
            field=models.DateField(null=True, default=None),
        ),
        migrations.AlterField(
            model_name='userlibpaper',
            name='date_last_modified',
            field=models.DateField(null=True, default=None),
        ),
    ]
