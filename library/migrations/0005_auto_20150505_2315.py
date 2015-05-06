# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_auto_20150505_2235'),
    ]

    operations = [
        migrations.AlterField(
            model_name='creator',
            name='email',
            field=models.EmailField(blank=True, null=True, max_length=254, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='e_issn',
            field=models.CharField(blank=True, max_length=10, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='ext_id',
            field=models.CharField(blank=True, max_length=30, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='issn',
            field=models.CharField(blank=True, max_length=10, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='language',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(null=True, default='', to='library.Publisher'),
        ),
        migrations.AlterField(
            model_name='journal',
            name='short_title',
            field=models.CharField(blank=True, max_length=100, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='title',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='url',
            field=models.URLField(blank=True, null=True, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='abstract',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='doi',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='ext_id',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='issue',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='journal',
            field=models.ForeignKey(null=True, default='', blank=True, to='library.Journal'),
        ),
        migrations.AlterField(
            model_name='paper',
            name='page',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='title',
            field=models.CharField(blank=True, max_length=500, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='url',
            field=models.URLField(blank=True, null=True, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='volume',
            field=models.CharField(blank=True, max_length=200, default=''),
        ),
    ]
