# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0002_auto_20160712_0926'),
    ]

    operations = [
        migrations.AddField(
            model_name='popover',
            name='priority',
            field=models.PositiveIntegerField(default=1),
        ),
    ]
