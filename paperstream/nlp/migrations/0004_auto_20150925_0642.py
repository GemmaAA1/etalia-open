# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import paperstream.nlp.mixins


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0003_auto_20150923_1734'),
    ]

    operations = [
        migrations.CreateModel(
            name='MostSimilar',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('model', models.OneToOneField(to='nlp.Model', related_name='ms')),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, paperstream.nlp.mixins.S3Mixin),
        ),
        migrations.RemoveField(
            model_name='papervectors',
            name='is_in_full_lsh',
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='model',
            field=models.ForeignKey(to='nlp.Model', default=7),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='time_lapse',
            field=models.IntegerField(default=-1, verbose_name='Days from right now', choices=[(7, '1 Week'), (30, '1 Month'), (60, '2 Months'), (365, '1 Year'), (-1, 'All')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('time_lapse', 'paper', 'model')]),
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='lsh',
        ),
    ]
