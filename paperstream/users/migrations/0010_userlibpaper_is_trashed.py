# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0009_userlibpaper_paper_provider_id'),
    ]

    operations = [
        migrations.AddField(
            model_name='userlibpaper',
            name='is_trashed',
            field=models.BooleanField(default=False),
        ),
    ]
