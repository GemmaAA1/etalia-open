# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0048_auto_20150516_0101'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', max_length=20, blank=True, choices=[('ZOTE', 'Zotero'), ('', 'Unknown'), ('IEEE', 'IEEE'), ('ELSV', 'Elsevier'), ('MEND', 'Mendeley'), ('ARXI', 'Arxiv'), ('PUMD', 'PubMed')]),
        ),
    ]
