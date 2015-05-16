# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0043_auto_20150516_0026'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', max_length=20, blank=True, choices=[('ELSV', 'Elsevier'), ('MEND', 'Mendeley'), ('PUMD', 'PubMed'), ('IEEE', 'IEEE'), ('ARXI', 'Arxiv'), ('ZOTE', 'Zotero'), ('', 'Unknown')]),
        ),
    ]
