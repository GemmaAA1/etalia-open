# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0010_auto_20150709_1953'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='model',
            options={'ordering': ['name']},
        ),
    ]
