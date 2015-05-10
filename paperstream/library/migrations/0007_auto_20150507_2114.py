# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150507_2011'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(null=True, max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(null=True, max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(null=True, max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(null=True, max_length=32, blank=True, db_index=True),
        ),
    ]
