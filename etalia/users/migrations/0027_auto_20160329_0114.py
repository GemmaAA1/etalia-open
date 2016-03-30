# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0026_auto_20160325_2053'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userthread',
            name='thread',
        ),
        migrations.RemoveField(
            model_name='userthread',
            name='user',
        ),
        migrations.RemoveField(
            model_name='user',
            name='threads',
        ),
        migrations.DeleteModel(
            name='UserThread',
        ),
    ]
