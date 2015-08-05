# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0028_auto_20150805_0704'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lsh',
            name='model',
            field=models.ForeignKey(to='nlp.Model', related_name='lsh'),
        ),
    ]
