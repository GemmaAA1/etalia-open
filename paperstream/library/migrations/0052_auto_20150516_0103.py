# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0051_auto_20150516_0103'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', max_length=20, blank=True, choices=[('PUMD', 'PubMed'), ('MEND', 'Mendeley'), ('ARXI', 'Arxiv'), ('ELSV', 'Elsevier'), ('', 'Unknown'), ('IEEE', 'IEEE'), ('ZOTE', 'Zotero')]),
        ),
    ]
