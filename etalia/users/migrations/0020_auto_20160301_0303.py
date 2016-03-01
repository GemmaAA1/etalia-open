# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0019_auto_20160212_0219'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='librarylayout',
            name='user',
        ),
        migrations.RemoveField(
            model_name='streamlayout',
            name='user',
        ),
        migrations.RemoveField(
            model_name='trendlayout',
            name='user',
        ),
        migrations.DeleteModel(
            name='LibraryLayout',
        ),
        migrations.DeleteModel(
            name='StreamLayout',
        ),
        migrations.DeleteModel(
            name='TrendLayout',
        ),
    ]
