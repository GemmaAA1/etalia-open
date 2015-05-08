# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0005_auto_20150507_0721'),
    ]

    operations = [
        migrations.AddField(
            model_name='journal',
            name='id_arx',
            field=models.CharField(max_length=32, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AddField(
            model_name='journal',
            name='id_issn',
            field=models.CharField(max_length=9, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AddField(
            model_name='journal',
            name='id_oth',
            field=models.CharField(max_length=9, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(max_length=32, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(max_length=32, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(max_length=32, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(max_length=32, db_index=True, blank=True, null=True, unique=True),
        ),
        migrations.AlterUniqueTogether(
            name='journal',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='journal',
            name='id_key',
        ),
        migrations.RemoveField(
            model_name='journal',
            name='id_val',
        ),
    ]
