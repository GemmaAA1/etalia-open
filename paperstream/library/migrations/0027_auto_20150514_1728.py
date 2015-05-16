# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import model_utils.fields
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0026_auto_20150514_1657'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='paper',
            name='is_aip',
        ),
        migrations.RemoveField(
            model_name='paper',
            name='is_pre_print',
        ),
        migrations.AddField(
            model_name='paper',
            name='publish_status',
            field=models.CharField(blank=True, default='', choices=[('ppublished', 'Paper Print'), ('epublished', 'e-Print only'), ('aheadofprint', 'e-Print ahead'), ('preprint', 'pre-Print'), ('', 'Unknown')], max_length=20),
        ),
        migrations.AddField(
            model_name='paper',
            name='publish_status_changed',
            field=model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='publish_status'),
        ),
        migrations.AddField(
            model_name='paper',
            name='source_changed',
            field=model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='source'),
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, default='', choices=[('ZOTE', 'Zotero'), ('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE'), ('MEND', 'Mendeley'), ('', 'Unknown'), ('PUMD', 'PubMed')], max_length=20),
        ),
    ]
