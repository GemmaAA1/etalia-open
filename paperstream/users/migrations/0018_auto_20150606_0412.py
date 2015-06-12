# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0017_remove_userlib_feed_status'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userlib',
            old_name='library_status',
            new_name='status',
        ),
    ]
