# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0012_auto_20160526_2358'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperengine',
            name='model',
            field=models.ForeignKey(related_name='ms', to='nlp.Model'),
        ),
    ]
