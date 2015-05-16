# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0036_auto_20150515_2256'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('MEND', 'Mendeley'), ('IEEE', 'IEEE'), ('PUMD', 'PubMed'), ('ARXI', 'Arxiv'), ('ELSV', 'Elsevier'), ('ZOTE', 'Zotero'), ('', 'Unknown')], blank=True, max_length=20, default=''),
        ),
    ]
