# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0013_auto_20160801_2322'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userlibpaper',
            name='date_created',
            field=models.DateField(db_index=True, null=True, default=None),
        ),
    ]
