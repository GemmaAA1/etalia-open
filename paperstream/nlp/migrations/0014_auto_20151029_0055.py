# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0013_mostsimilar_journal_ratio'),
    ]

    operations = [
        migrations.AddField(
            model_name='mostsimilar',
            name='is_active',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='mostsimilar',
            name='model',
            field=models.ForeignKey(related_name='ms', to='nlp.Model'),
        ),
        migrations.AlterUniqueTogether(
            name='mostsimilar',
            unique_together=set([('model', 'journal_ratio')]),
        ),
    ]
