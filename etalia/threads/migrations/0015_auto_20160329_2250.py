# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0014_auto_20160329_2248'),
    ]

    operations = [
        migrations.RenameField(
            model_name='threadcomment',
            old_name='author',
            new_name='user',
        ),
        migrations.RenameField(
            model_name='threadpost',
            old_name='author',
            new_name='user',
        ),
        migrations.AlterUniqueTogether(
            name='threadcomment',
            unique_together=set([('post', 'user', 'content')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadpost',
            unique_together=set([('thread', 'content', 'user')]),
        ),
    ]
