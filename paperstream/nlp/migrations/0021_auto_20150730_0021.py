# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0020_auto_20150729_1818'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='paperneighbors',
            name='model',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='neighbor1',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='neighbor2',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='neighbor3',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='neighbor4',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='neighbor5',
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='lsh',
            field=models.ForeignKey(to='nlp.LSH', default=3),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='neighbors',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), null=True, size=10),
        ),
    ]
