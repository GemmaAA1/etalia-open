# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0002_auto_20160602_0013'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='threadfeed',
            name='thread_scores',
        ),
        migrations.RemoveField(
            model_name='threadfeed',
            name='user',
        ),
        migrations.AlterUniqueTogether(
            name='threadscore',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='threadscore',
            name='thread',
        ),
        migrations.RemoveField(
            model_name='threadscore',
            name='thread_feed',
        ),
        migrations.DeleteModel(
            name='ThreadFeed',
        ),
        migrations.DeleteModel(
            name='ThreadScore',
        ),
    ]
