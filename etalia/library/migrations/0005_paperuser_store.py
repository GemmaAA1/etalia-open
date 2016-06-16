# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_paperuserhistory'),
    ]

    operations = [
        migrations.AddField(
            model_name='paperuser',
            name='store',
            field=models.PositiveIntegerField(choices=[(1, 'Pinned'), (2, 'Trashed')], default=None, null=True),
        ),
    ]
