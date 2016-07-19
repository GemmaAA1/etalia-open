# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0006_auto_20160714_1534'),
    ]

    operations = [
        migrations.AlterField(
            model_name='popover',
            name='template_path',
            field=models.CharField(blank=True, max_length=128, null=True),
        ),
        migrations.AlterField(
            model_name='popover',
            name='title',
            field=models.CharField(blank=True, max_length=256, null=True),
        ),
    ]
