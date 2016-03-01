# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0016_auto_20160113_0003'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='stream_author_weight',
            field=models.FloatField(verbose_name='Author weight', default=1.0),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_journal_weight',
            field=models.FloatField(verbose_name='Journal weight', default=1.0),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_vector_weight',
            field=models.FloatField(verbose_name='Content weight', default=1.0),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_altmetric_weight',
            field=models.FloatField(verbose_name='Altmetric weight', default=1.0),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_doc_weight',
            field=models.FloatField(verbose_name='Content weight', default=1.0),
        ),
    ]
