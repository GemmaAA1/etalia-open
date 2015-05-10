# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='author',
            name='email',
            field=models.EmailField(max_length=254, default='', blank=True),
        ),
        migrations.AlterField(
            model_name='authorposition',
            name='position',
            field=models.IntegerField(default=1),
        ),
        migrations.AlterField(
            model_name='journal',
            name='url',
            field=models.URLField(default='', blank=True),
        ),
        migrations.AlterField(
            model_name='paper',
            name='url',
            field=models.URLField(default='', blank=True),
        ),
    ]
