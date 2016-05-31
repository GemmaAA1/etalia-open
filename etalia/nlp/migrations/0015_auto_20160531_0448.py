# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0014_paperengine_journal_boost'),
    ]

    operations = [
        migrations.AlterField(
            model_name='threadengine',
            name='model',
            field=models.ForeignKey(related_name='te', to='nlp.Model'),
        ),
    ]
