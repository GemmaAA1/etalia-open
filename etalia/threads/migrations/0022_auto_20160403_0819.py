# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0021_auto_20160403_0817'),
    ]

    operations = [
        migrations.AlterField(
            model_name='threaduserhistory',
            name='date',
            field=models.DateTimeField(auto_now_add=True),
        ),
    ]
