# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0005_auto_20160712_1020'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='popover',
            name='body',
        ),
        migrations.AddField(
            model_name='popover',
            name='template_path',
            field=models.CharField(default='', max_length=128),
            preserve_default=False,
        ),
    ]
