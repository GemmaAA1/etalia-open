# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0013_auto_20160329_0742'),
    ]

    operations = [
        migrations.RenameField(
            model_name='thread',
            old_name='owner',
            new_name='user',
        ),
        migrations.AlterUniqueTogether(
            name='thread',
            unique_together=set([('type', 'user', 'title', 'paper')]),
        ),
        migrations.RemoveField(
            model_name='thread',
            name='users',
        ),
    ]
