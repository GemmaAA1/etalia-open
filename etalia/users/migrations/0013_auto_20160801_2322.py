# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0012_userlib_state'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='user',
            name='is_alpha',
        ),
        migrations.AddField(
            model_name='user',
            name='type',
            field=models.IntegerField(choices=[(1, 'Individual'), (2, 'Third-Party')], verbose_name='Type', default=1),
        ),
    ]
