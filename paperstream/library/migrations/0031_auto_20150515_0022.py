# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0030_auto_20150515_0022'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', max_length=20, blank=True, choices=[('ARXI', 'Arxiv'), ('PUMD', 'PubMed'), ('ELSV', 'Elsevier'), ('ZOTE', 'Zotero'), ('IEEE', 'IEEE'), ('MEND', 'Mendeley'), ('', 'Unknown')]),
        ),
    ]
