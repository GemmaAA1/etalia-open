# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0002_auto_20151103_1907'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mostsimilar',
            name='journal_ratio',
            field=models.FloatField(choices=[('0.0', '0 %'), ('0.05', '5 %'), ('0.10', '10 %'), ('0.15', '15 %'), ('0.20', '20 %'), ('0.25', '25 %'), ('0.30', '30 %')], default=0.0),
        ),
    ]
