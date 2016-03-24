# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0001_initial'),
    ]

    operations = [
        migrations.RenameField(
            model_name='thread',
            old_name='author',
            new_name='user',
        ),
        migrations.RenameField(
            model_name='threadpost',
            old_name='author',
            new_name='user',
        ),
        migrations.RenameField(
            model_name='threadpostcomment',
            old_name='author',
            new_name='user',
        ),
        migrations.AlterField(
            model_name='thread',
            name='content',
            field=models.TextField(blank=True, default='', verbose_name='Content', null=True),
        ),
        migrations.AlterField(
            model_name='thread',
            name='paper',
            field=models.ForeignKey(blank=True, null=True, to='library.Paper', default=None, verbose_name='Related Paper'),
        ),
        migrations.AlterField(
            model_name='thread',
            name='title',
            field=models.CharField(max_length=256, verbose_name='Title'),
        ),
        migrations.AlterField(
            model_name='thread',
            name='type',
            field=models.IntegerField(default=1, choices=[(1, 'Question'), (2, 'Paper')], verbose_name='Type'),
        ),
    ]
