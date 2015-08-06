# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0029_auto_20150805_1748'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([('journal', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([('paper', 'model')]),
        ),
    ]
