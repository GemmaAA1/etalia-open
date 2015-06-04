# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0015_userstats'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userlibjournal',
            name='score',
        ),
    ]
