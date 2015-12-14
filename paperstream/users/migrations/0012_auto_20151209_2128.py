# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0011_auto_20151209_2112'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userlibjournal',
            old_name='papers_in_journal',
            new_name='occurrence',
        ),
        migrations.AlterModelOptions(
            name='userlibjournal',
            options={'ordering': ('-occurrence',)},
        ),
    ]
