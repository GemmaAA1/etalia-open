# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0002_auto_20150930_0451'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='usertaste',
            unique_together=set([('user', 'paper')]),
        ),
    ]
