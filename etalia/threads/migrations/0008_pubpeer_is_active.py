# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0007_auto_20160831_1923'),
    ]

    operations = [
        migrations.AddField(
            model_name='pubpeer',
            name='is_active',
            field=models.BooleanField(default=True),
        ),
    ]
