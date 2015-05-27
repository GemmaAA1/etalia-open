# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='consumerjournalstat',
            name='date',
        ),
        migrations.AddField(
            model_name='consumerjournalstat',
            name='datetime',
            field=models.DateTimeField(default=datetime.datetime(2015, 5, 21, 23, 14, 13, 428862, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='consumerjournalstat',
            name='number_papers_fetched',
            field=models.IntegerField(),
        ),
        migrations.AlterField(
            model_name='consumerjournalstat',
            name='number_papers_recorded',
            field=models.IntegerField(),
        ),
    ]
