# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0006_consumerjournal_date_last_cons'),
    ]

    operations = [
        migrations.RenameField(
            model_name='consumer',
            old_name='max_ret',
            new_name='ret_max',
        ),
    ]
