# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0013_auto_20150511_1836'),
    ]

    operations = [
        migrations.RenameField(
            model_name='paper',
            old_name='is_locked',
            new_name='is_trusted',
        ),
    ]
