# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0008_auto_20150522_2007'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userlibjournal',
            old_name='user',
            new_name='userlib',
        ),
        migrations.RenameField(
            model_name='userlibpaper',
            old_name='user',
            new_name='userlib',
        ),
        migrations.RenameField(
            model_name='userlibstats',
            old_name='user',
            new_name='userlib',
        ),
        migrations.AlterField(
            model_name='userlib',
            name='user',
            field=models.OneToOneField(to=settings.AUTH_USER_MODEL, related_name='userlib', serialize=False, primary_key=True),
        ),
    ]
