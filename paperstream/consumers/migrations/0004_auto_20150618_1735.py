# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0003_auto_20150618_0625'),
    ]

    operations = [
        migrations.RenameField(
            model_name='consumerjournal',
            old_name='base_countdown_day',
            new_name='base_coundown_period',
        ),
        migrations.RenameField(
            model_name='consumerjournal',
            old_name='countdown_day',
            new_name='coundown_period',
        ),
    ]
