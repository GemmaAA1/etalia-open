# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0007_auto_20160720_2223'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usersettings',
            name='trend_altmetric_weight',
            field=models.FloatField(verbose_name='Altmetric weight', default=1),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_doc_weight',
            field=models.FloatField(verbose_name='Title/Abstract weight', default=0.2),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_score_threshold',
            field=models.FloatField(verbose_name='Specificity', default=0.1),
        ),
    ]
