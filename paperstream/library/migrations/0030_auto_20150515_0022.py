# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0029_auto_20150514_2142'),
    ]

    operations = [
        migrations.RenameField(
            model_name='paper',
            old_name='terms',
            new_name='key_terms',
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, max_length=20, choices=[('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('ZOTE', 'Zotero'), ('MEND', 'Mendeley'), ('', 'Unknown'), ('PUMD', 'PubMed'), ('IEEE', 'IEEE')], default=''),
        ),
    ]
