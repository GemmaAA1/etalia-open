# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0026_auto_20160421_0046'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='thread',
            options={'ordering': ('-published_at', '-modified')},
        ),
    ]
