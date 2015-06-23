# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0003_auto_20150618_1612'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userfeed',
            old_name='paper_out',
            new_name='paper_match',
        ),
        migrations.RenameField(
            model_name='userfeed',
            old_name='paper_in',
            new_name='paper_seed',
        ),
    ]
