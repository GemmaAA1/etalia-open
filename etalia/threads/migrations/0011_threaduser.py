# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('threads', '0010_auto_20160325_1830'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadUser',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('is_pinned', models.BooleanField(default=False)),
                ('is_banned', models.BooleanField(default=False)),
                ('is_joined', models.BooleanField(default=False)),
                ('is_left', models.BooleanField(default=False)),
                ('first_joined_at', models.DateTimeField(null=True, blank=True, default=None)),
                ('last_left_at', models.DateTimeField(null=True, blank=True, default=None)),
                ('num_comments', models.PositiveIntegerField(default=0)),
                ('thread', models.ForeignKey(to='threads.Thread')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-num_comments', 'first_joined_at'],
            },
        ),
    ]
