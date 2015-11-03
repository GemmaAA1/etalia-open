# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0002_auto_20151102_0359'),
    ]

    operations = [
        migrations.AlterField(
            model_name='trend',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='trends'),
        ),
    ]
