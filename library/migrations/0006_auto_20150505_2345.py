# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0005_auto_20150505_2315'),
    ]

    operations = [
        migrations.RenameField(
            model_name='journal',
            old_name='paper_counter',
            new_name='lib_size',
        ),
    ]
