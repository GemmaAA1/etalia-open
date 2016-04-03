# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0020_auto_20160403_0743'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='threaduserhistory',
            name='action',
        ),
        migrations.AddField(
            model_name='threaduserhistory',
            name='difference',
            field=models.CharField(max_length=256, default=''),
        ),
    ]
