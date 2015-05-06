# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0009_auto_20150506_0212'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='journal',
            options={'ordering': ['title']},
        ),
    ]
