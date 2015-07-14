# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_paper_id_isbn'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='date_ep',
            field=models.DateField(blank=True, null=True, default=None, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='date_pp',
            field=models.DateField(blank=True, null=True, default=None, db_index=True),
        ),
    ]
