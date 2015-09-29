# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_paper_date_fs'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='date_fs',
            field=models.DateField(auto_now_add=True, db_index=True),
        ),
    ]
