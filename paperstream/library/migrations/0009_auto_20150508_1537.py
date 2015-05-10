# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0008_auto_20150508_1533'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(db_index=True, unique=True, null=True, max_length=32, blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(db_index=True, unique=True, null=True, max_length=32, blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(db_index=True, unique=True, null=True, max_length=32, blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(db_index=True, unique=True, null=True, max_length=32, blank=True),
        ),
    ]
