# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0022_auto_20150731_2301'),
    ]

    operations = [
        migrations.AddField(
            model_name='model',
            name='fields',
            field=models.ManyToManyField(related_name='fields', to='nlp.FieldUseInModel'),
        ),
        migrations.AlterUniqueTogether(
            name='fielduseinmodel',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='fielduseinmodel',
            name='model',
        ),
    ]
