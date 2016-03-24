# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0005_thread_privacy'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='threadpost',
            unique_together=set([('thread', 'content', 'user')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadpostcomment',
            unique_together=set([('post', 'user', 'content')]),
        ),
        migrations.RemoveField(
            model_name='threadpost',
            name='position',
        ),
        migrations.RemoveField(
            model_name='threadpostcomment',
            name='position',
        ),
    ]
