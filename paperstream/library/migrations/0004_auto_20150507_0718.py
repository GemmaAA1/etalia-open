# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20150507_0644'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='paper',
            name='identifiers',
        ),
        migrations.AddField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(max_length=32, unique=True, null=True),
        ),
        migrations.AddField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(max_length=32, unique=True, null=True),
        ),
        migrations.AddField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(max_length=32, unique=True, null=True),
        ),
        migrations.AddField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(max_length=32, unique=True, null=True),
        ),
    ]
