# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0011_auto_20160526_2315'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paperengine',
            name='model',
            field=models.ForeignKey(unique=True, to='nlp.Model', related_name='pe'),
        ),
        migrations.AlterUniqueTogether(
            name='paperengine',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='paperengine',
            name='journal_ratio',
        ),
    ]
