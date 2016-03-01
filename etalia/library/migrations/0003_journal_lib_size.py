# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_remove_journal_lib_size'),
    ]

    operations = [
        migrations.AddField(
            model_name='journal',
            name='lib_size',
            field=models.IntegerField(default=0),
        ),
    ]
