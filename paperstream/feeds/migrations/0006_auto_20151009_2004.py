# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0005_auto_20151009_1916'),
    ]

    operations = [
        migrations.AlterField(
            model_name='discoverfeed',
            name='user',
            field=models.OneToOneField(to=settings.AUTH_USER_MODEL, related_name='discover'),
        ),
    ]
