# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0003_auto_20151103_0101'),
    ]

    operations = [
        migrations.RenameField(
            model_name='trendmatches',
            old_name='trend_feed',
            new_name='trend',
        ),
        migrations.AlterUniqueTogether(
            name='trendmatches',
            unique_together=set([('trend', 'paper')]),
        ),
    ]
