# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0002_auto_20160602_2323'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='usertaste',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='usertaste',
            name='paper',
        ),
        migrations.RemoveField(
            model_name='usertaste',
            name='user',
        ),
        migrations.DeleteModel(
            name='UserTaste',
        ),
    ]
