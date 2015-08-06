# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150710_0053'),
        ('nlp', '0018_auto_20150723_1810'),
    ]

    operations = [
        migrations.CreateModel(
            name='PaperNeighbors',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='lsh',
            name='time_lapse',
            field=models.IntegerField(default=None, verbose_name='In the past for', null=True, blank=True, choices=[(61, '2 Months')]),
        ),
        migrations.AlterField(
            model_name='model',
            name='hs',
            field=models.IntegerField(default=1),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='neighbor1',
            field=models.ForeignKey(to='library.Paper', related_name='neighbor1'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='neighbor2',
            field=models.ForeignKey(to='library.Paper', related_name='neighbor2'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='neighbor3',
            field=models.ForeignKey(to='library.Paper', related_name='neighbor3'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='neighbor4',
            field=models.ForeignKey(to='library.Paper', related_name='neighbor4'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='neighbor5',
            field=models.ForeignKey(to='library.Paper', related_name='neighbor5'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
    ]
