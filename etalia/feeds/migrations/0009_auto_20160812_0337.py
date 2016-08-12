# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0008_auto_20160812_0332'),
    ]

    operations = [
        migrations.AlterField(
            model_name='streampapers',
            name='score',
            field=models.FloatField(default=0.0, db_index=True),
        ),
        migrations.AlterField(
            model_name='threadfeedthreads',
            name='score',
            field=models.FloatField(default=0.0, db_index=True),
        ),
        migrations.AlterField(
            model_name='trendpapers',
            name='score',
            field=models.FloatField(default=0.0, db_index=True),
        ),
    ]
