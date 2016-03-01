# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0017_auto_20160120_0109'),
    ]

    operations = [
        migrations.RenameField(
            model_name='usertaste',
            old_name='is_ticked',
            new_name='is_banned',
        ),
        migrations.RenameField(
            model_name='usertaste',
            old_name='is_liked',
            new_name='is_pinned',
        ),
    ]
