# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0006_auto_20161013_1712'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='PubPeerConsumer',
            new_name='ConsumerPubPeer',
        ),
    ]
