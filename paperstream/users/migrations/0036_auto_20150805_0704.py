# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0035_auto_20150805_0538'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='userlibjournal',
            unique_together=set([('userlib', 'journal')]),
        ),
    ]
