# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0054_auto_20150516_0103'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', choices=[('ELSV', 'Elsevier'), ('ZOTE', 'Zotero'), ('', 'Unknown'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE'), ('PUMD', 'PubMed'), ('MEND', 'Mendeley')], max_length=20, blank=True),
        ),
    ]
