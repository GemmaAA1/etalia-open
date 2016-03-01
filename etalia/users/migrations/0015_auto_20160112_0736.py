# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0014_auto_20151218_0821'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='trend_method_args',
            field=jsonfield.fields.JSONField(null=True, default=None, blank=True),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_method',
            field=models.IntegerField(verbose_name='Method', choices=[(0, 'Content Based Scoring (simple)')], default=0),
        ),
    ]
