# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0034_auto_20150515_1620'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', choices=[('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('MEND', 'Mendeley'), ('ZOTE', 'Zotero'), ('PUMD', 'PubMed'), ('IEEE', 'IEEE'), ('', 'Unknown')], max_length=20, blank=True),
        ),
    ]
