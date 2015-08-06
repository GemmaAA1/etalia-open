# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0026_auto_20150802_0643'),
    ]

    operations = [
        migrations.AlterField(
            model_name='lsh',
            name='state',
            field=models.CharField(max_length=3, default='None', choices=[('NON', 'None'), ('BUS', 'Busy'), ('IDL', 'Idle')]),
        ),
        migrations.AlterUniqueTogether(
            name='lsh',
            unique_together=set([('model', 'time_lapse')]),
        ),
    ]
