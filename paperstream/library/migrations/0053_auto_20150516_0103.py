# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0052_auto_20150516_0103'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', blank=True, max_length=20, choices=[('IEEE', 'IEEE'), ('PUMD', 'PubMed'), ('', 'Unknown'), ('MEND', 'Mendeley'), ('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('ZOTE', 'Zotero')]),
        ),
    ]
