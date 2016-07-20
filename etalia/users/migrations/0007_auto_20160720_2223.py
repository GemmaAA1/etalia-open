# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0006_auto_20160719_0124'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='stream_score_threshold',
            field=models.FloatField(verbose_name='Specificity', default=0.2),
        ),
        migrations.AddField(
            model_name='usersettings',
            name='trend_score_threshold',
            field=models.FloatField(verbose_name='Specificity', default=0.2),
        ),
        migrations.AlterField(
            model_name='user',
            name='init_step',
            field=models.CharField(max_length=3, default='NON', help_text='Tag where init user stands', choices=[('NON', 'uninitialized'), ('LIB', 'Syncing library'), ('STR', 'Syncing feed Papers'), ('TRE', 'Syncing feed Trend'), ('THR', 'Syncing feed Thread'), ('POP', 'Syncing popovers'), ('IDL', 'done')]),
        ),
    ]
