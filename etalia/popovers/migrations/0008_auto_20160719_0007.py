# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0007_auto_20160718_2241'),
    ]

    operations = [
        migrations.AddField(
            model_name='userpopover',
            name='display',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='popover',
            name='priority',
            field=models.PositiveIntegerField(validators=[django.core.validators.MinValueValidator(1), django.core.validators.MaxValueValidator(9)], default=1),
        ),
    ]
