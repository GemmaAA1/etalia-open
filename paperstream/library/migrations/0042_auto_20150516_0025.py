# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0041_auto_20150516_0025'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', blank=True, choices=[('ZOTE', 'Zotero'), ('IEEE', 'IEEE'), ('ARXI', 'Arxiv'), ('ELSV', 'Elsevier'), ('PUMD', 'PubMed'), ('MEND', 'Mendeley'), ('', 'Unknown')], max_length=20),
        ),
    ]
