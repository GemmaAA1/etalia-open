# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0023_auto_20150513_2113'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('IEEE', 'IEEE'), ('PBMD', 'PubMed'), ('ELSV', 'Elsevier'), ('', 'Unknown'), ('ZOTE', 'Zotero'), ('MEND', 'Mendeley'), ('ARXI', 'Arxiv')], max_length=4, blank=True),
        ),
    ]
