# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0021_auto_20170604_2140'),
    ]

    operations = [
        migrations.AddField(
            model_name='userlibpaper',
            name='is_orcid',
            field=models.BooleanField(default=False),
        ),
    ]
