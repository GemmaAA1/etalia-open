# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0011_auto_20150511_1735'),
    ]

    operations = [
        migrations.RenameField(
            model_name='journal',
            old_name='is_valid',
            new_name='is_locked',
        ),
    ]
