# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0014_auto_20150814_2102'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userfeed',
            name='user',
            field=models.ForeignKey(related_name='feeds', to=settings.AUTH_USER_MODEL),
        ),
    ]
