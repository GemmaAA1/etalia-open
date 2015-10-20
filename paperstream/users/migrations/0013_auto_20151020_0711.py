# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0012_userfeedlayout'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userfeedlayout',
            old_name='stream',
            new_name='stream_filter',
        ),
        migrations.RenameField(
            model_name='userfeedlayout',
            old_name='trend',
            new_name='trend_filter',
        ),
        migrations.AlterField(
            model_name='userfeedlayout',
            name='user',
            field=models.OneToOneField(to=settings.AUTH_USER_MODEL),
        ),
    ]
