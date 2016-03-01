# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0018_auto_20160211_0014'),
    ]

    operations = [
        migrations.RenameField(
            model_name='usertaste',
            old_name='context_source',
            new_name='source',
        ),
    ]
