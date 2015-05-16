# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0028_auto_20150514_2139'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='paper',
            options={},
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, max_length=20, choices=[('IEEE', 'IEEE'), ('MEND', 'Mendeley'), ('ZOTE', 'Zotero'), ('ELSV', 'Elsevier'), ('', 'Unknown'), ('ARXI', 'Arxiv'), ('PUMD', 'PubMed')], default=''),
        ),
    ]
