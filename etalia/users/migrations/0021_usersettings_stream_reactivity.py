# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0020_auto_20160301_0303'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='stream_reactivity',
            field=models.FloatField(default=1.0, verbose_name='Reactivity'),
        ),
    ]
