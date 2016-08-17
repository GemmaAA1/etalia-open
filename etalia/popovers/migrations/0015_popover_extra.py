# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0014_auto_20160816_0523'),
    ]

    operations = [
        migrations.AddField(
            model_name='popover',
            name='extra',
            field=models.CharField(max_length=64, blank=True, default=''),
        ),
    ]
