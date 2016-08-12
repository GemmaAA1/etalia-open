# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0007_auto_20160721_1745'),
    ]

    operations = [
        migrations.AlterField(
            model_name='streampapers',
            name='date',
            field=models.DateField(db_index=True),
        ),
        migrations.AlterField(
            model_name='threadfeedthreads',
            name='date',
            field=models.DateField(db_index=True),
        ),
        migrations.AlterField(
            model_name='trendpapers',
            name='date',
            field=models.DateField(db_index=True),
        ),
    ]
