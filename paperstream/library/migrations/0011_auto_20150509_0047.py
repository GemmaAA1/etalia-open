# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0010_auto_20150508_1541'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='id_arx',
            field=models.CharField(default='', max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_eissn',
            field=models.CharField(default='', max_length=9, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_issn',
            field=models.CharField(default='', max_length=9, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='id_oth',
            field=models.CharField(default='', max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_arx',
            field=models.CharField(default='', max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_doi',
            field=models.CharField(default='', max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_oth',
            field=models.CharField(default='', max_length=32, blank=True, db_index=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pmi',
            field=models.CharField(default='', max_length=32, blank=True, db_index=True),
        ),
    ]
