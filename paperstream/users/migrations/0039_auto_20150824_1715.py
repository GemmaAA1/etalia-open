# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0038_auto_20150814_2102'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userlibpaper',
            name='paper',
            field=models.ForeignKey(related_name='userlib_paper', to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='userlibpaper',
            name='userlib',
            field=models.ForeignKey(related_name='userlib_paper', to='users.UserLib'),
        ),
    ]
