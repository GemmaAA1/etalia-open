# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0029_remove_user_photo'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='usersettings',
            name='stream_method_args',
        ),
        migrations.RemoveField(
            model_name='usersettings',
            name='stream_narrowness',
        ),
        migrations.RemoveField(
            model_name='usersettings',
            name='stream_time_lapse',
        ),
        migrations.RemoveField(
            model_name='usersettings',
            name='trend_method_args',
        ),
        migrations.RemoveField(
            model_name='usersettings',
            name='trend_narrowness',
        ),
        migrations.RemoveField(
            model_name='usersettings',
            name='trend_time_lapse',
        ),
    ]
