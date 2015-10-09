# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_stats'),
    ]

    operations = [
        migrations.AddField(
            model_name='stats',
            name='nb_papers_last_month',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='stats',
            name='nb_papers_last_two_week',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='stats',
            name='nb_papers_last_week',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='stats',
            name='nb_papers_last_year',
            field=models.IntegerField(default=0),
        ),
    ]
