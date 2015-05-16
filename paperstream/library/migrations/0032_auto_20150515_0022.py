# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0031_auto_20150515_0022'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', choices=[('ELSV', 'Elsevier'), ('PUMD', 'PubMed'), ('ZOTE', 'Zotero'), ('ARXI', 'Arxiv'), ('MEND', 'Mendeley'), ('IEEE', 'IEEE'), ('', 'Unknown')], blank=True, max_length=20),
        ),
    ]
