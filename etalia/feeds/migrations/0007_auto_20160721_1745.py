# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0006_auto_20160720_2223'),
    ]

    operations = [
        migrations.RenameField(
            model_name='stream',
            old_name='last_update',
            new_name='updated_at',
        ),
        migrations.RenameField(
            model_name='threadfeed',
            old_name='last_update',
            new_name='updated_at',
        ),
        migrations.RenameField(
            model_name='trend',
            old_name='last_update',
            new_name='updated_at',
        ),
    ]
