# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0009_auto_20150716_2107'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userfeed',
            name='papers_match',
            field=models.ManyToManyField(related_name='papers_match', through='feeds.UserFeedPaper', to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='userfeed',
            name='papers_seed',
            field=models.ManyToManyField(related_name='papers_seed', to='library.Paper'),
        ),
        migrations.AlterField(
            model_name='userfeedvector',
            name='feed',
            field=models.ForeignKey(related_name='vectors', to='feeds.UserFeed'),
        ),
    ]
