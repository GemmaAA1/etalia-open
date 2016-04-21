# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0025_auto_20160420_1838'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='threadneighbor',
            name='thread',
        ),
        migrations.DeleteModel(
            name='ThreadNeighbor',
        ),
    ]
