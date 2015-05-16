# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0046_auto_20150516_0053'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(default='', choices=[('PUMD', 'PubMed'), ('ELSV', 'Elsevier'), ('ZOTE', 'Zotero'), ('', 'Unknown'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE'), ('MEND', 'Mendeley')], blank=True, max_length=20),
        ),
    ]
