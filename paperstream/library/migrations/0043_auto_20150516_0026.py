# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0042_auto_20150516_0025'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', max_length=20, blank=True, choices=[('IEEE', 'IEEE'), ('ZOTE', 'Zotero'), ('ARXI', 'Arxiv'), ('', 'Unknown'), ('ELSV', 'Elsevier'), ('PUMD', 'PubMed'), ('MEND', 'Mendeley')]),
        ),
    ]
