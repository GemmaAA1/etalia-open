# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0007_auto_20150513_1830'),
    ]

    operations = [
        migrations.AddField(
            model_name='consumer',
            name='type',
            field=models.CharField(default='PBMD', choices=[('PBMD', 'PubMed'), ('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE')], max_length=4),
            preserve_default=False,
        ),
    ]
