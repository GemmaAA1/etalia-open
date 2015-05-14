# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0021_auto_20150513_2110'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('MEND', 'Mendeley'), ('', 'Unknown'), ('IEEE', 'IEEE'), ('ZOTE', 'Zotero'), ('ELSV', 'Elsevier'), ('PBMD', 'PubMed'), ('ARXI', 'Arxiv')], blank=True, max_length=4),
        ),
    ]
