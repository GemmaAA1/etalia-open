# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0006_auto_20150708_2248'),
    ]

    operations = [
        migrations.RenameField(
            model_name='papervectors',
            old_name='mod',
            new_name='model',
        ),
    ]
