# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0016_auto_20151027_0521'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userfeedlayout',
            name='library_filter',
            field=jsonfield.fields.JSONField(null=True),
        ),
        migrations.AlterField(
            model_name='userfeedlayout',
            name='stream_filter',
            field=jsonfield.fields.JSONField(null=True),
        ),
        migrations.AlterField(
            model_name='userfeedlayout',
            name='trend_filter',
            field=jsonfield.fields.JSONField(null=True),
        ),
    ]
