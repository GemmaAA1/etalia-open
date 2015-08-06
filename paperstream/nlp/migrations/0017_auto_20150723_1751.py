# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0016_auto_20150723_1750'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lshfmodel',
            name='model',
            field=models.OneToOneField(to='nlp.Model', related_name='lshf'),
        ),
    ]
