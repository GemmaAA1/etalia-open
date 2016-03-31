# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0002_auto_20160323_0531'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='thread',
            unique_together=set([('type', 'user', 'title', 'paper')]),
        ),
    ]
