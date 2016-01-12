# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0006_streammatches_date'),
    ]

    operations = [
        migrations.AddField(
            model_name='trendmatches',
            name='date',
            field=models.DateField(default=datetime.datetime(2016, 1, 7, 16, 1, 11, 866630, tzinfo=utc)),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='streammatches',
            name='date',
            field=models.DateField(default=datetime.date(2016, 1, 7)),
            preserve_default=False,
        ),
    ]
