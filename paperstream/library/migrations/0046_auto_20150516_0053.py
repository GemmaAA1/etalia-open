# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0045_auto_20150516_0050'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('ELSV', 'Elsevier'), ('IEEE', 'IEEE'), ('MEND', 'Mendeley'), ('ARXI', 'Arxiv'), ('PUMD', 'PubMed'), ('ZOTE', 'Zotero'), ('', 'Unknown')], max_length=20, blank=True, default=''),
        ),
    ]
