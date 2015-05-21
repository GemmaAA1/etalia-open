# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0002_auto_20150520_0558'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumerjournal',
            name='last_date_cons',
            field=models.DateTimeField(null=True, default=datetime.datetime(2015, 3, 22, 6, 27, 38, 515534, tzinfo=utc)),
        ),
    ]
