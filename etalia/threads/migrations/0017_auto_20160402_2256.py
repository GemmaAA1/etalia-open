# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0016_auto_20160401_0007'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='threaduser',
            unique_together=set([('thread', 'user')]),
        ),
    ]
