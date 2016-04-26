# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0027_auto_20160421_1748'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='thread',
            unique_together=set([]),
        ),
    ]
