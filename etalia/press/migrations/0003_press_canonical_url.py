# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('press', '0002_auto_20160908_1806'),
    ]

    operations = [
        migrations.AddField(
            model_name='press',
            name='canonical_url',
            field=models.URLField(null=True, blank=True, default=''),
        ),
    ]
