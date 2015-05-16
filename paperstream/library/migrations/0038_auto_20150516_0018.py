# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0037_auto_20150516_0017'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', blank=True, choices=[('', 'Unknown'), ('MEND', 'Mendeley'), ('IEEE', 'IEEE'), ('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('PUMD', 'PubMed'), ('ZOTE', 'Zotero')], max_length=20),
        ),
    ]
