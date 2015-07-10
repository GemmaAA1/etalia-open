# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0026_auto_20150709_2119'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='model',
            field=models.ForeignKey(default=23, to='nlp.Model'),
            preserve_default=False,
        ),
    ]
