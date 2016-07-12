# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0003_popover_priority'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='userpopover',
            unique_together=set([('popover', 'status')]),
        ),
    ]
