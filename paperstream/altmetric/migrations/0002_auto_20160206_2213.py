# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('altmetric', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='altmetricmodel',
            name='image',
            field=models.URLField(default=''),
        ),
        migrations.AddField(
            model_name='altmetricmodel',
            name='type',
            field=models.CharField(default='', max_length=32),
        ),
    ]
