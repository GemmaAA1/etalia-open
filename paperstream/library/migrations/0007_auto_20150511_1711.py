# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150511_0838'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='period',
            field=models.CharField(choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular'), (None, 'Unknown')], max_length=200, blank=True, null=True, default=None),
        ),
    ]
