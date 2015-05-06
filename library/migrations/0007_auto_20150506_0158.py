# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150505_2345'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='journal',
            field=models.ForeignKey(null=True, to='library.Journal'),
        ),
    ]
