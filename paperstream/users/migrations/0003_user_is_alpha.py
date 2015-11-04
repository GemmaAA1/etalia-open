# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0002_auto_20151103_1907'),
    ]

    operations = [
        migrations.AddField(
            model_name='user',
            name='is_alpha',
            field=models.BooleanField(verbose_name='alpha', help_text='Designates whether this user should be treated as an early adopter user.', default=True),
        ),
    ]
