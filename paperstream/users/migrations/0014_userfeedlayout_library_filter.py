# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0013_auto_20151020_0711'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfeedlayout',
            name='library_filter',
            field=jsonfield.fields.JSONField(default={}),
            preserve_default=False,
        ),
    ]
