# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_auto_20150507_0718'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(default=None, max_length=32, unique=True, null=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(default=None, max_length=32, unique=True, null=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(default=None, max_length=32, unique=True, null=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(default=None, max_length=32, unique=True, null=True),
        ),
    ]
