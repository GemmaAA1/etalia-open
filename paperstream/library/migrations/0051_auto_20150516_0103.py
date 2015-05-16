# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0050_auto_20150516_0101'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(max_length=20, default='', blank=True, choices=[('IEEE', 'IEEE'), ('ELSV', 'Elsevier'), ('MEND', 'Mendeley'), ('ZOTE', 'Zotero'), ('', 'Unknown'), ('ARXI', 'Arxiv'), ('PUMD', 'PubMed')]),
        ),
    ]
