# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0030_auto_20150715_1731'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='method',
            field=models.IntegerField(default=1, choices=[(1, 'Average'), (2, 'Average threshold'), (3, 'Average date weighted')], verbose_name='Scoring Algo'),
        ),
    ]
