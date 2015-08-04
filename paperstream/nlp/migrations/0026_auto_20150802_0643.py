# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0025_auto_20150731_2316'),
    ]

    operations = [
        migrations.AlterField(
            model_name='model',
            name='_size',
            field=models.IntegerField(validators=[django.core.validators.MaxValueValidator(300)], default=128),
        ),
    ]
