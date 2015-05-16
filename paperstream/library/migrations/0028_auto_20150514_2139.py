# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0027_auto_20150514_1728'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='paper',
            options={'ordering': ['-created']},
        ),
        migrations.RenameField(
            model_name='paper',
            old_name='date_p',
            new_name='date_pp',
        ),
        migrations.AlterField(
            model_name='paper',
            name='publish_status',
            field=models.CharField(blank=True, default='', choices=[('ppublish', 'Paper Print'), ('epublish', 'e-Print only'), ('aheadofprint', 'e-Print ahead'), ('preprint', 'pre-Print'), ('', 'Unknown')], max_length=20),
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, default='', choices=[('', 'Unknown'), ('ARXI', 'Arxiv'), ('MEND', 'Mendeley'), ('PUMD', 'PubMed'), ('ELSV', 'Elsevier'), ('ZOTE', 'Zotero'), ('IEEE', 'IEEE')], max_length=20),
        ),
    ]
