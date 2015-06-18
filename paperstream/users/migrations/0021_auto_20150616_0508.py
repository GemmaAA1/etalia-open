# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import users.validators


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0020_auto_20150609_1817'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='first_name',
            field=models.CharField(blank=True, default='', max_length=255, validators=[users.validators.validate_first_name]),
        ),
        migrations.AlterField(
            model_name='user',
            name='last_name',
            field=models.CharField(blank=True, default='', max_length=255, validators=[users.validators.validate_last_name]),
        ),
    ]
