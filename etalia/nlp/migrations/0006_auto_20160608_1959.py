# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0005_auto_20160608_1954'),
    ]

    operations = [
        migrations.RenameField(
            model_name='threadengine',
            old_name='paper_boost',
            new_name='score_paper_boost',
        ),
        migrations.RenameField(
            model_name='threadengine',
            old_name='user_boost',
            new_name='score_user_boost',
        ),
    ]
