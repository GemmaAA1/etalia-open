# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0037_auto_20150806_0308'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='scoring_method',
            field=models.IntegerField(choices=[(1, 'Average'), (2, 'Average Threshold binary'), (3, 'Average weighted journal'), (4, 'Average weighted journal-date')], default=1, verbose_name='Scoring Algo'),
        ),
        migrations.AlterUniqueTogether(
            name='userlibpaper',
            unique_together=set([('userlib', 'paper')]),
        ),
    ]
