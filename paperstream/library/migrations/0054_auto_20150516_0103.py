# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0053_auto_20150516_0103'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, default='', choices=[('', 'Unknown'), ('MEND', 'Mendeley'), ('ELSV', 'Elsevier'), ('IEEE', 'IEEE'), ('ARXI', 'Arxiv'), ('PUMD', 'PubMed'), ('ZOTE', 'Zotero')], max_length=20),
        ),
    ]
