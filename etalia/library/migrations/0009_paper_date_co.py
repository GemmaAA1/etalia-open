# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0008_auto_20160812_0538'),
    ]

    operations = [
        migrations.AddField(
            model_name='paper',
            name='date_co',
            field=models.DateField(null=True, blank=True, db_index=True),
        ),
    ]
