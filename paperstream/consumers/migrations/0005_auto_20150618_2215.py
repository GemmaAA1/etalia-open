# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0004_auto_20150618_1735'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='consumerjournal',
            options={'ordering': ['last_date_cons']},
        ),
        migrations.AlterUniqueTogether(
            name='consumerjournal',
            unique_together=set([]),
        ),
    ]
