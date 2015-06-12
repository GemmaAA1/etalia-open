# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0016_remove_userlibjournal_score'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userlib',
            name='feed_status',
        ),
    ]
