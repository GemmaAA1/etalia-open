# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mostsimilar',
            name='journal_ratio',
            field=models.FloatField(default=0.0, choices=[('0.0', '0 %'), ('0.1', '10 %'), ('0.15', '15 %'), ('0.20', '20 %'), ('0.25', '25 %'), ('0.30', '30 %')]),
        ),
    ]
