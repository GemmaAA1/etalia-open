# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0009_auto_20160422_2358'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='journalneighbors',
            name='journal',
        ),
        migrations.RemoveField(
            model_name='journalneighbors',
            name='pe',
        ),
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='journalvectors',
            name='journal',
        ),
        migrations.RemoveField(
            model_name='journalvectors',
            name='model',
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='pe',
        ),
        migrations.RemoveField(
            model_name='paperneighbors',
            name='paper',
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='papervectors',
            name='model',
        ),
        migrations.RemoveField(
            model_name='papervectors',
            name='paper',
        ),
        migrations.AlterUniqueTogether(
            name='threadneighbors',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='threadneighbors',
            name='pe',
        ),
        migrations.RemoveField(
            model_name='threadneighbors',
            name='thread',
        ),
        migrations.AlterUniqueTogether(
            name='threadvectors',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='threadvectors',
            name='model',
        ),
        migrations.RemoveField(
            model_name='threadvectors',
            name='thread',
        ),
        migrations.DeleteModel(
            name='JournalNeighbors',
        ),
        migrations.DeleteModel(
            name='JournalVectors',
        ),
        migrations.DeleteModel(
            name='PaperNeighbors',
        ),
        migrations.DeleteModel(
            name='PaperVectors',
        ),
        migrations.DeleteModel(
            name='ThreadNeighbors',
        ),
        migrations.DeleteModel(
            name='ThreadVectors',
        ),
    ]
