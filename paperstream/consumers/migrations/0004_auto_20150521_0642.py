# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0003_auto_20150521_0627'),
    ]

    operations = [
        migrations.RenameField(
            model_name='consumerjournal',
            old_name='last_number_papers_retrieved',
            new_name='last_number_papers_recorded',
        ),
        migrations.AlterField(
            model_name='consumerjournal',
            name='last_date_cons',
            field=models.DateTimeField(null=True, default=datetime.datetime(2015, 3, 22, 6, 42, 51, 69791, tzinfo=utc)),
        ),
    ]
