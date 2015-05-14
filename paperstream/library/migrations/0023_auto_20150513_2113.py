# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0022_auto_20150513_2112'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('IEEE', 'IEEE'), ('', 'Unknown'), ('ZOTE', 'Zotero'), ('ARXI', 'Arxiv'), ('MEND', 'Mendeley'), ('ELSV', 'Elsevier'), ('PBMD', 'PubMed')], max_length=4, blank=True),
        ),
    ]
