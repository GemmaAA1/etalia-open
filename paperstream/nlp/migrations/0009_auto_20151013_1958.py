# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0008_auto_20151013_1935'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mostsimilarstatus',
            name='ms',
            field=models.ForeignKey(related_name='status', to='nlp.MostSimilar'),
        ),
    ]
