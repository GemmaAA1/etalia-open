# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0025_usersettings'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='model',
            field=models.ForeignKey(null=True, to='nlp.Model'),
        ),
    ]
