# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0002_auto_20150513_0524'),
    ]

    operations = [
        migrations.AddField(
            model_name='consumer',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 5, 13, 5, 27, 4, 742956, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='consumer',
            name='modified',
            field=models.DateTimeField(auto_now=True, default=datetime.datetime(2015, 5, 13, 5, 27, 11, 846960, tzinfo=utc)),
            preserve_default=False,
        ),
    ]
