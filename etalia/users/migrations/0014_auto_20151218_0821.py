# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0013_auto_20151218_0158'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='stream_method_args',
            field=jsonfield.fields.JSONField(default=None, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_method',
            field=models.IntegerField(verbose_name='Method', choices=[(0, 'Content Based Simple')], default=0),
        ),
    ]
