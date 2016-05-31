# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0003_auto_20151104_2202'),
    ]

    operations = [
        migrations.AddField(
            model_name='journalneighbors',
            name='pe',
            field=models.ForeignKey(to='nlp.models.PaperEngine', default=5),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='pe',
            field=models.ForeignKey(to='nlp.models.PaperEngine', default=5),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='paperneighbors',
            name='time_lapse',
            field=models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year'), (-1, 'All')], verbose_name='Days from right now'),
        ),
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([('journal', 'pe')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('time_lapse', 'paper', 'pe')]),
        ),
        migrations.RemoveField(
            model_name='journalneighbors',
            name='model',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='model',
        ),
    ]
