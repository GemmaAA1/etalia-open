# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_auto_20150505_2202'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='e_issn',
            field=models.CharField(null=True, blank=True, max_length=10),
        ),
        migrations.AlterField(
            model_name='journal',
            name='ext_id',
            field=models.CharField(null=True, blank=True, max_length=30),
        ),
        migrations.AlterField(
            model_name='journal',
            name='issn',
            field=models.CharField(null=True, blank=True, max_length=10),
        ),
        migrations.AlterField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(null=True, to='library.Publisher', blank=True),
        ),
        migrations.AlterField(
            model_name='journal',
            name='title',
            field=models.CharField(null=True, blank=True, max_length=200),
        ),
    ]
