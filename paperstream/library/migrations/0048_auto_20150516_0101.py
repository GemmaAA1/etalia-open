# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0047_auto_20150516_0100'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, choices=[('IEEE', 'IEEE'), ('PUMD', 'PubMed'), ('ELSV', 'Elsevier'), ('MEND', 'Mendeley'), ('', 'Unknown'), ('ARXI', 'Arxiv'), ('ZOTE', 'Zotero')], max_length=20, default=''),
        ),
    ]
