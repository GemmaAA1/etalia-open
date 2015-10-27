# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0014_userfeedlayout_library_filter'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='usertaste',
            unique_together=set([('user', 'paper', 'scoring_method')]),
        ),
    ]
