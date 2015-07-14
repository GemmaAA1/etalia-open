# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0006_auto_20150623_1656'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userfeedpaper',
            name='is_in_user_lib',
        ),
        migrations.AddField(
            model_name='userfeedpaper',
            name='is_score_computed',
            field=models.BooleanField(default=False),
        ),
    ]
