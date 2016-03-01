# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('altmetric', '0002_auto_20160206_2213'),
    ]

    operations = [
        migrations.AlterField(
            model_name='altmetricmodel',
            name='type',
            field=models.CharField(default='zzzzzzzz', max_length=8),
        ),
    ]
