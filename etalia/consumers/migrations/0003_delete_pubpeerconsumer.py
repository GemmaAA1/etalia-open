# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0002_auto_20160823_2343'),
    ]

    operations = [
        migrations.DeleteModel(
            name='PubPeerConsumer',
        ),
    ]
