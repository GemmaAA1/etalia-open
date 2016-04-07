# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0018_auto_20160402_2354'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='threaduser',
            options={'ordering': ['-num_comments', 'first_joined_at']},
        ),
    ]
