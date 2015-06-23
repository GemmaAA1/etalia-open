# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import feeds.validators


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0005_auto_20150623_1549'),
    ]

    operations = [
        migrations.RenameField(
            model_name='userfeed',
            old_name='paper_match',
            new_name='papers_match',
        ),
        migrations.RenameField(
            model_name='userfeed',
            old_name='paper_seed',
            new_name='papers_seed',
        ),
        migrations.AlterField(
            model_name='userfeed',
            name='name',
            field=models.CharField(default='main', max_length=100, validators=[feeds.validators.validate_feed_name]),
        ),
    ]
