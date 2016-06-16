# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0003_auto_20160603_1807'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='userlibauthor',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='userlibauthor',
            name='author',
        ),
        migrations.RemoveField(
            model_name='userlibauthor',
            name='userlib',
        ),
        migrations.AlterUniqueTogether(
            name='userlibjournal',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='userlibjournal',
            name='journal',
        ),
        migrations.RemoveField(
            model_name='userlibjournal',
            name='userlib',
        ),
        migrations.RemoveField(
            model_name='userlib',
            name='authors',
        ),
        migrations.RemoveField(
            model_name='userlib',
            name='journals',
        ),
        migrations.RemoveField(
            model_name='userlibpaper',
            name='is_trashed',
        ),
        migrations.DeleteModel(
            name='UserLibAuthor',
        ),
        migrations.DeleteModel(
            name='UserLibJournal',
        ),
    ]
