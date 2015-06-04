# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0013_auto_20150528_0459'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userstats',
            name='user',
        ),
        migrations.DeleteModel(
            name='UserStats',
        ),
    ]
