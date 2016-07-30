# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0009_usersession'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='usersession',
            name='session',
        ),
        migrations.RemoveField(
            model_name='usersession',
            name='user',
        ),
        migrations.DeleteModel(
            name='UserSession',
        ),
    ]
