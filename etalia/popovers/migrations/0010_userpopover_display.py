# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0009_remove_userpopover_display'),
    ]

    operations = [
        migrations.AddField(
            model_name='userpopover',
            name='display',
            field=models.PositiveIntegerField(choices=[(2, 'Display'), (1, 'Hide')], default=1),
        ),
    ]
