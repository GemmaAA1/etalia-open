# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0003_auto_20160603_0740'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='threaduserhistory',
            options={'ordering': ('-created',)},
        ),
    ]
