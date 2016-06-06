# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0005_paperuser_store'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='paperuserhistory',
            options={'ordering': ('-created',)},
        ),
    ]
