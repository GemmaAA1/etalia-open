# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0005_auto_20150511_0218'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='period',
            field=models.CharField(default='IRR', choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular')], blank=True, null=True, max_length=200),
        ),
    ]