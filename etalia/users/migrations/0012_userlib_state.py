# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0011_remove_userlib_state'),
    ]

    operations = [
        migrations.AddField(
            model_name='userlib',
            name='state',
            field=models.IntegerField(choices=[(1, 'Uninitialized'), (2, 'Idle'), (3, 'Syncing'), (4, 'Need New OAuth')], default=1),
        ),
    ]
