# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0007_auto_20160422_2319'),
    ]

    operations = [
        migrations.AddField(
            model_name='threadneighbors',
            name='ms',
            field=models.ForeignKey(default=1, to='nlp.models.ThreadEngine'),
            preserve_default=False,
        ),
    ]
