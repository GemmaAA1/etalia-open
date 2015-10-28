# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0005_auto_20150924_1741'),
    ]

    operations = [
        migrations.AddField(
            model_name='consumerjournal',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 10, 27, 23, 0, 24, 375488, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='consumerjournal',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 10, 27, 23, 0, 28, 359780, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='consumerjournalstat',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 10, 27, 23, 0, 34, 631749, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='consumerjournalstat',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 10, 27, 23, 0, 38, 319898, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
    ]
