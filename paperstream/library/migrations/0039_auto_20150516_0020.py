# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0038_auto_20150516_0018'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(choices=[('MEND', 'Mendeley'), ('ZOTE', 'Zotero'), ('ARXI', 'Arxiv'), ('ELSV', 'Elsevier'), ('IEEE', 'IEEE'), ('PUMD', 'PubMed'), ('', 'Unknown')], default='', blank=True, max_length=20),
        ),
    ]
