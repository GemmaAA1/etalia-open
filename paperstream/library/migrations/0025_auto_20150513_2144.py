# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0024_auto_20150513_2144'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(max_length=4, choices=[('ZOTE', 'Zotero'), ('IEEE', 'IEEE'), ('ELSV', 'Elsevier'), ('MEND', 'Mendeley'), ('PBMD', 'PubMed'), ('ARXI', 'Arxiv'), ('', 'Unknown')], blank=True),
        ),
    ]
