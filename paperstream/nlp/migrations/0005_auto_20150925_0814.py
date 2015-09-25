# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0004_auto_20150925_0642'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='lsh',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='lsh',
            name='model',
        ),
        migrations.DeleteModel(
            name='LSH',
        ),
    ]
