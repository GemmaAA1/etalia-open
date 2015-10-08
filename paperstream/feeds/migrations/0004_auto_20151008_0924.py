# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0003_auto_20150930_0553'),
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
