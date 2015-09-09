# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0001_initial'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='consumerjournal',
            unique_together=set([('journal', 'consumer')]),
        ),
    ]
