# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0007_auto_20160323_0856'),
    ]

    operations = [
        migrations.RenameField(
            model_name='thread',
            old_name='user',
            new_name='owner',
        ),
        migrations.RenameField(
            model_name='threadmember',
            old_name='user',
            new_name='member',
        ),
        migrations.RenameField(
            model_name='threadpost',
            old_name='user',
            new_name='author',
        ),
        migrations.RenameField(
            model_name='threadpostcomment',
            old_name='user',
            new_name='author',
        ),
        migrations.AlterUniqueTogether(
            name='thread',
            unique_together=set([('type', 'owner', 'title', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadpost',
            unique_together=set([('thread', 'content', 'author')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadpostcomment',
            unique_together=set([('post', 'author', 'content')]),
        ),
    ]
