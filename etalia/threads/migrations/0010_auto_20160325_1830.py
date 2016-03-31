# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0009_auto_20160325_0642'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='threadmember',
            name='member',
        ),
        migrations.RemoveField(
            model_name='threadmember',
            name='thread',
        ),
        migrations.RemoveField(
            model_name='threaduserstate',
            name='thread',
        ),
        migrations.RemoveField(
            model_name='threaduserstate',
            name='user',
        ),
        migrations.RemoveField(
            model_name='thread',
            name='members',
        ),
        migrations.AlterField(
            model_name='thread',
            name='owner',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='threads_owned'),
        ),
        migrations.DeleteModel(
            name='ThreadMember',
        ),
        migrations.DeleteModel(
            name='ThreadUserState',
        ),
    ]
