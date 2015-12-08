# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0004_auto_20151105_0422'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='stream_narrowness',
            field=models.IntegerField(choices=[(-2, 'Narrowest'), (-1, 'Narrower'), (0, 'Normal'), (1, 'Broader'), (2, 'Broadest')], default=0, verbose_name='Narrowness'),
        ),
        migrations.AddField(
            model_name='usersettings',
            name='trend_narrowness',
            field=models.IntegerField(choices=[(-2, 'Narrowest'), (-1, 'Narrower'), (0, 'Normal'), (1, 'Broader'), (2, 'Broadest')], default=0, verbose_name='Narrowness'),
        ),
    ]
