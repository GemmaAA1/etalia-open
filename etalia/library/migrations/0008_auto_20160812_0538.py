# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0007_auto_20160606_0555'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperuser',
            name='store',
            field=models.PositiveIntegerField(choices=[(1, 'Added'), (2, 'Trashed')], null=True, db_index=True, default=None),
        ),
        migrations.AlterField(
            model_name='paperuser',
            name='watch',
            field=models.PositiveIntegerField(choices=[(1, 'Pinned'), (2, 'Banned')], null=True, db_index=True, default=None),
        ),
    ]
