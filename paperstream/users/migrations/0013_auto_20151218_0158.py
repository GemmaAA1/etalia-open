# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0012_auto_20151209_2128'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_method',
            field=models.IntegerField(default=0, verbose_name='Method', choices=[(0, 'Content Based Simple'), (1, 'Max'), (2, 'Average'), (3, 'Average Threshold binary'), (4, 'Average weighted journal'), (5, 'Average weighted journal-date'), (6, 'Occurrence Count')]),
        ),
    ]
