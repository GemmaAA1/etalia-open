# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0005_auto_20150521_2307'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userstats',
            old_name='doing',
            new_name='state',
        ),
        migrations.RemoveField(
            model_name='userstats',
            name='when',
        ),
        migrations.AddField(
            model_name='userstats',
            name='datetime',
            field=models.DateTimeField(default=datetime.datetime(2015, 5, 21, 23, 14, 31, 371451, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
    ]
