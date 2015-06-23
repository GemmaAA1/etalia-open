# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0002_auto_20150616_0531'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfeed',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 11, 45, 285565, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userfeed',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 11, 51, 53886, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userfeedpaper',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 6, 44385, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userfeedpaper',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 17, 97990, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
    ]
