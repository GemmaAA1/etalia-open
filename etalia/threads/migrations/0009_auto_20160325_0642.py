# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0008_auto_20160323_1853'),
    ]

    operations = [
        migrations.RenameField(
            model_name='threaduserstate',
            old_name='is_added',
            new_name='is_joined',
        ),
        migrations.RenameField(
            model_name='threaduserstate',
            old_name='is_trashed',
            new_name='is_left',
        ),
    ]
