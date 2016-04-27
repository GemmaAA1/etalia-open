# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0027_auto_20160329_0114'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='relationship',
            unique_together=set([('from_user', 'to_user')]),
        ),
    ]
