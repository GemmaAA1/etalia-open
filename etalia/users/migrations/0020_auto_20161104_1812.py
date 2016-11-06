# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0019_auto_20160816_0703'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='threadfeed_score_threshold',
            field=models.FloatField(verbose_name='Specificity', default=-1),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_score_threshold',
            field=models.FloatField(verbose_name='Specificity', default=-1),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_score_threshold',
            field=models.FloatField(verbose_name='Specificity', default=-1),
        ),
    ]
