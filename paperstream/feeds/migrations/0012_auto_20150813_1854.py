# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0011_auto_20150716_2122'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='userfeed',
            unique_together=set([('name', 'user')]),
        ),
    ]
