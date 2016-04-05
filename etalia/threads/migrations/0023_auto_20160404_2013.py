# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0022_auto_20160403_0819'),
    ]

    operations = [
        migrations.AlterField(
            model_name='thread',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, null=True, related_name='threads_owned'),
        ),
    ]
