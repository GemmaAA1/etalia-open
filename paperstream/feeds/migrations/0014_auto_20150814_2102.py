# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0013_auto_20150813_1937'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='userfeedpaper',
            unique_together=set([('feed', 'paper')]),
        ),
    ]
