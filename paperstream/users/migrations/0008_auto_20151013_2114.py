# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0007_auto_20151009_2115'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='userlibjournal',
            options={'ordering': ('-papers_in_journal',)},
        ),
    ]
