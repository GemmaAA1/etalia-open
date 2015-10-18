# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0008_auto_20151013_2114'),
    ]

    operations = [
        migrations.AddField(
            model_name='userlibpaper',
            name='paper_provider_id',
            field=models.CharField(default='', max_length=64),
        ),
    ]
