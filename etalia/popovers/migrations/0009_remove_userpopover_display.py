# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0008_auto_20160719_0007'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userpopover',
            name='display',
        ),
    ]
