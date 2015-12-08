# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0005_auto_20151207_0234'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='email_digest_frequency',
            field=models.IntegerField(default=7, choices=[(7, 'Weekly'), (15, 'Bi-Weekly'), (30, 'Monthly'), (30, 'Bi-Monthly'), (-1, 'Never')], verbose_name='Email digest frequency'),
        ),
    ]
