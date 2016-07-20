# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0011_userpopoverupdatedisplay'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='userpopover',
            options={'ordering': ('-popover__type', 'popover__priority')},
        ),
    ]
