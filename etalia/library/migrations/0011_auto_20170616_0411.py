# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0010_journal_is_in_fixture'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='date_fs',
            field=models.DateField(db_index=True, null=True, blank=True),
        ),
    ]
