# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0020_auto_20150513_2048'),
    ]

    operations = [
        migrations.AddField(
            model_name='paper',
            name='terms',
            field=models.TextField(default='', blank=True),
        ),
        migrations.AddField(
            model_name='paper',
            name='type',
            field=models.CharField(max_length=4, choices=[('JOUR', 'Journal Article'), ('LETT', 'Letter'), ('EDIT', 'Editorial'), ('NEWS', 'News'), ('PROC', 'Proceedings'), ('REVI', 'Review'), ('DRAF', 'Draft'), ('', 'Unknown')], blank=True, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(max_length=4, choices=[('PBMD', 'PubMed'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE'), ('ELSV', 'Elsevier'), ('', 'Unknown'), ('MEND', 'Mendeley'), ('ZOTE', 'Zotero')], blank=True),
        ),
    ]
