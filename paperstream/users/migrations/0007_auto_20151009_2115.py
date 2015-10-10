# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0006_auto_20151008_0924'),
    ]

    operations = [
        migrations.RenameField(
            model_name='usertaste',
            old_name='is_disliked',
            new_name='is_ticked',
        ),
    ]
