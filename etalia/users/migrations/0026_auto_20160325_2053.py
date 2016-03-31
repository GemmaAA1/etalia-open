# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0025_auto_20160325_1830'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='userthread',
            options={'ordering': ['-num_comments', 'first_joined_at']},
        ),
        migrations.RenameField(
            model_name='userthread',
            old_name='joined_at',
            new_name='first_joined_at',
        ),
        migrations.RenameField(
            model_name='userthread',
            old_name='left_at',
            new_name='last_left_at',
        ),
    ]
