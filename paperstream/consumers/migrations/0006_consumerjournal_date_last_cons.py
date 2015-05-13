# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0005_auto_20150513_0744'),
    ]

    operations = [
        migrations.AddField(
            model_name='consumerjournal',
            name='date_last_cons',
            field=models.DateTimeField(null=True, default=None),
        ),
    ]
