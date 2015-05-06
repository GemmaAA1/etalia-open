# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20150505_2232'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='e_issn',
            field=models.CharField(null=True, max_length=10, default='', blank=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='ext_id',
            field=models.CharField(null=True, max_length=30, default='', blank=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='issn',
            field=models.CharField(null=True, max_length=10, default='', blank=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='title',
            field=models.CharField(null=True, max_length=200, default='', blank=True),
        ),
    ]
