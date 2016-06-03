# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_roll_back_deltatime',
            field=models.IntegerField(null=True, default=None, verbose_name='Roll-back time (months)', blank=True),
        ),
    ]
