# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0017_auto_20150819_1840'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userfeed',
            name='papers_seed',
        ),
    ]
