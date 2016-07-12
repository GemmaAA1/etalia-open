# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='popover',
            name='anchor',
            field=models.CharField(null=True, blank=True, max_length=128),
        ),
    ]
