# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20160605_2237'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperuser',
            name='store',
            field=models.PositiveIntegerField(null=True, choices=[(1, 'Added'), (2, 'Trashed')], default=None),
        ),
    ]
