# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0023_auto_20150731_2309'),
    ]

    operations = [
        migrations.AlterField(
            model_name='model',
            name='fields',
            field=models.ManyToManyField(to='nlp.FieldUseInModel', null=True, related_name='fields'),
        ),
    ]
