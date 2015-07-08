# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import users.validators


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0023_auto_20150618_1612'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='email',
            field=models.EmailField(unique=True, max_length=255, db_index=True, verbose_name='Email'),
        ),
        migrations.AlterField(
            model_name='user',
            name='first_name',
            field=models.CharField(max_length=255, default='', verbose_name='First Name', blank=True, validators=[users.validators.validate_first_name]),
        ),
        migrations.AlterField(
            model_name='user',
            name='last_name',
            field=models.CharField(max_length=255, default='', verbose_name='Last Name', blank=True, validators=[users.validators.validate_last_name]),
        ),
        migrations.AlterField(
            model_name='user',
            name='username',
            field=models.CharField(max_length=255, db_index=True, verbose_name='username (UNUSED)', blank=True, default=''),
        ),
    ]
