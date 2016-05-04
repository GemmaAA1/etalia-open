# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0023_auto_20160308_0835'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_author_weight',
            field=models.FloatField(default=0.5, verbose_name='Author weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_journal_weight',
            field=models.FloatField(default=0.5, verbose_name='Journal weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_roll_back_deltatime',
            field=models.IntegerField(default=36, verbose_name='Roll-back time (months)'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_vector_weight',
            field=models.FloatField(default=0.5, verbose_name='Content weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_altmetric_weight',
            field=models.FloatField(default=0.5, verbose_name='Altmetric weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_doc_weight',
            field=models.FloatField(default=0.5, verbose_name='Content weight'),
        ),
    ]
