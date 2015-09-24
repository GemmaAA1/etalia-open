# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0004_auto_20150924_1728'),
    ]

    operations = [
        migrations.AlterField(
            model_name='consumerjournal',
            name='last_date_cons',
            field=models.DateTimeField(blank=True, default=None, null=True),
        ),
    ]
