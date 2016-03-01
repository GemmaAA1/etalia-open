# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0007_auto_20160107_1601'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='trendmatches',
            options={'ordering': ['-score']},
        ),
    ]
