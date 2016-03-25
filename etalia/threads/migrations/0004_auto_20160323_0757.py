# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0003_auto_20160323_0552'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='threadmember',
            options={'ordering': ['-num_comments', 'first_joined_at']},
        ),
    ]
