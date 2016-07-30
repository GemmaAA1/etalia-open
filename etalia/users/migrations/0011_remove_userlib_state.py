# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0010_auto_20160729_2157'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userlib',
            name='state',
        ),
    ]
