# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='affiliation',
            name='city',
            field=models.CharField(default='', blank=True, max_length=50),
        ),
        migrations.AlterField(
            model_name='affiliation',
            name='country',
            field=models.CharField(default='', blank=True, max_length=50),
        ),
        migrations.AlterField(
            model_name='affiliation',
            name='department',
            field=models.CharField(default='', blank=True, max_length=200),
        ),
        migrations.AlterField(
            model_name='affiliation',
            name='institution',
            field=models.CharField(default='', blank=True, max_length=200),
        ),
        migrations.AlterField(
            model_name='affiliation',
            name='state',
            field=models.CharField(default='', blank=True, max_length=20),
        ),
    ]
