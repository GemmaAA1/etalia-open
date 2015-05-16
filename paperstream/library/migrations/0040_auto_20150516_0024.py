# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0039_auto_20150516_0020'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(max_length=20, blank=True, default='', choices=[('IEEE', 'IEEE'), ('ARXI', 'Arxiv'), ('', 'Unknown'), ('PUMD', 'PubMed'), ('MEND', 'Mendeley'), ('ZOTE', 'Zotero'), ('ELSV', 'Elsevier')]),
        ),
    ]
