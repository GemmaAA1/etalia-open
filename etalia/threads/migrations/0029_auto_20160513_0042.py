# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0028_auto_20160425_1850'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='threaduserinvite',
            unique_together=set([('thread', 'from_user', 'to_user')]),
        ),
    ]
