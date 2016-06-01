# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0030_auto_20160526_0031'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='usersettings',
            name='stream_model',
        ),
        migrations.RemoveField(
            model_name='usersettings',
            name='trend_model',
        ),
    ]
