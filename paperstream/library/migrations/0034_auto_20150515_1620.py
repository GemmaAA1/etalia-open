# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import library.validators


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0033_auto_20150515_1554'),
    ]

    operations = [
        migrations.AlterField(
            model_name='author',
            name='last_name',
            field=models.CharField(validators=[library.validators.validate_author_names], max_length=100, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, choices=[('ARXI', 'Arxiv'), ('PUMD', 'PubMed'), ('IEEE', 'IEEE'), ('ZOTE', 'Zotero'), ('', 'Unknown'), ('MEND', 'Mendeley'), ('ELSV', 'Elsevier')], max_length=20, default=''),
        ),
    ]
