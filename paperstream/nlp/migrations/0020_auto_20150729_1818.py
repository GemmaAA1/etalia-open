# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0019_auto_20150729_1817'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperneighbors',
            name='neighbor1',
            field=models.ForeignKey(related_name='neighbor1', null=True, to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='neighbor2',
            field=models.ForeignKey(related_name='neighbor2', null=True, to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='neighbor3',
            field=models.ForeignKey(related_name='neighbor3', null=True, to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='neighbor4',
            field=models.ForeignKey(related_name='neighbor4', null=True, to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='neighbor5',
            field=models.ForeignKey(related_name='neighbor5', null=True, to='library.Paper'),
        ),
    ]
