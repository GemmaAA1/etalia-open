# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0007_auto_20150506_0158'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(default='', to='library.Publisher', null=True, blank=True),
        ),
    ]
