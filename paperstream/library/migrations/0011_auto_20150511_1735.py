# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0010_auto_20150511_1731'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(blank=True, null=True, to='library.Publisher', default=None),
        ),
    ]
