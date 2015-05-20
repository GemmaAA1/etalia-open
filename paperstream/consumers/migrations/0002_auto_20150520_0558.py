# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumerjournal',
            name='last_date_cons',
            field=models.DateTimeField(default=datetime.datetime(2015, 3, 21, 5, 58, 30, 327636, tzinfo=utc), null=True),
        ),
    ]
