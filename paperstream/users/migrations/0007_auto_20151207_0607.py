# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0006_usersettings_email_digest_frequency'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='trend_time_lapse',
            field=models.IntegerField(default=60, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year'), (-1, 'All')], verbose_name='Time range'),
        ),
    ]
