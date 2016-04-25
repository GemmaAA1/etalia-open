# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0008_threadneighbors_ms'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='threadneighbors',
            unique_together=set([('time_lapse', 'thread', 'ms')]),
        ),
    ]
