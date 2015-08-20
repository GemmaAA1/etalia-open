# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150710_0053'),
        ('feeds', '0018_remove_userfeed_papers_seed'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfeed',
            name='papers_seed',
            field=models.ManyToManyField(through='feeds.UserFeedSeedPaper', to='library.Paper', related_name='feed_seed'),
        ),
    ]
