# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0040_auto_20150516_0024'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('PUMD', 'PubMed'), ('ELSV', 'Elsevier'), ('IEEE', 'IEEE'), ('', 'Unknown'), ('MEND', 'Mendeley'), ('ARXI', 'Arxiv'), ('ZOTE', 'Zotero')], max_length=20, blank=True, default=''),
        ),
    ]
