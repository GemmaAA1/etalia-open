# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('threads', '0012_thread_users'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadComment',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('content', models.TextField(default='', blank=True)),
                ('author', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('post', models.ForeignKey(related_name='comments', to='threads.ThreadPost')),
            ],
            options={
                'ordering': ('created',),
            },
        ),
        migrations.AlterUniqueTogether(
            name='threadpostcomment',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='threadpostcomment',
            name='author',
        ),
        migrations.RemoveField(
            model_name='threadpostcomment',
            name='post',
        ),
        migrations.DeleteModel(
            name='ThreadPostComment',
        ),
        migrations.AlterUniqueTogether(
            name='threadcomment',
            unique_together=set([('post', 'author', 'content')]),
        ),
    ]
