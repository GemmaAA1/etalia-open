# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0010_auto_20160325_1830'),
        ('users', '0024_auto_20160322_2311'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserThread',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('is_pinned', models.BooleanField(default=False)),
                ('is_banned', models.BooleanField(default=False)),
                ('is_joined', models.BooleanField(default=False)),
                ('is_left', models.BooleanField(default=False)),
                ('joined_at', models.DateTimeField(blank=True, default=None, null=True)),
                ('left_at', models.DateTimeField(blank=True, default=None, null=True)),
                ('num_comments', models.PositiveIntegerField(default=0)),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-num_comments', 'joined_at'],
            },
        ),
        migrations.AddField(
            model_name='user',
            name='threads',
            field=models.ManyToManyField(to='threads.Thread', through='users.UserThread'),
        ),
    ]
