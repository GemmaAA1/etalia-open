# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0022_auto_20160307_2126'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='usersettings',
            name='stream_reactivity',
        ),
        migrations.AddField(
            model_name='usersettings',
            name='stream_roll_back_deltatime',
            field=models.IntegerField(default=36, verbose_name='Rolling back time (Lower time will focus on recent addition to your library'),
        ),
    ]
