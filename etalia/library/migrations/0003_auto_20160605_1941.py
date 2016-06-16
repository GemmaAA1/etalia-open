# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0002_auto_20160603_1807'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='paperuserhistory',
            name='paperuser',
        ),
        migrations.RemoveField(
            model_name='paperuserhistory',
            name='paperuser_ptr',
        ),
        migrations.DeleteModel(
            name='PaperUserHistory',
        ),
    ]
