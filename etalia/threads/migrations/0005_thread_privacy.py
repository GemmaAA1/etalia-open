# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0004_auto_20160323_0757'),
    ]

    operations = [
        migrations.AddField(
            model_name='thread',
            name='privacy',
            field=models.IntegerField(default=1, choices=[(1, 'Public'), (2, 'Private')]),
        ),
    ]
