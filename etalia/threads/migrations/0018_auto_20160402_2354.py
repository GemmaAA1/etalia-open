# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0017_auto_20160402_2256'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='threaduser',
            options={'permissions': (('change_threaduser', ''), ('delete_threaduser', '')), 'ordering': ['-num_comments', 'first_joined_at']},
        ),
    ]
