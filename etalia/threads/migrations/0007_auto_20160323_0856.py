# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0006_auto_20160323_0855'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='threadpost',
            options={'ordering': ('created',)},
        ),
        migrations.AlterModelOptions(
            name='threadpostcomment',
            options={'ordering': ('created',)},
        ),
    ]
