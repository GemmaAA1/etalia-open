# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0004_auto_20150619_0639'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfeedpaper',
            name='is_disliked',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='userfeedpaper',
            name='is_in_user_lib',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='userfeedpaper',
            name='is_liked',
            field=models.BooleanField(default=False),
        ),
    ]
