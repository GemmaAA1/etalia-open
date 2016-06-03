# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import etalia.core.mixins


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Thread',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('type', models.IntegerField(default=1, choices=[(1, 'Question'), (2, 'Paper')], verbose_name='Type')),
                ('privacy', models.IntegerField(default=1, choices=[(1, 'Public'), (2, 'Private')])),
                ('title', models.CharField(default='', verbose_name='Title', max_length=256)),
                ('content', models.TextField(blank=True, default='', verbose_name='Content', null=True)),
                ('published_at', models.DateTimeField(blank=True, null=True)),
            ],
            options={
                'ordering': ('-published_at', '-modified'),
            },
        ),
        migrations.CreateModel(
            name='ThreadComment',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('content', models.TextField(blank=True, default='')),
            ],
            options={
                'ordering': ('created',),
            },
        ),
        migrations.CreateModel(
            name='ThreadFeed',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='main', max_length=128)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('status', models.CharField(blank=True, default='NON', choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadPost',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('content', models.TextField(blank=True, default='')),
            ],
            options={
                'ordering': ('created',),
            },
        ),
        migrations.CreateModel(
            name='ThreadScore',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('new', models.BooleanField(default=True)),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.CreateModel(
            name='ThreadUser',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('watch', models.PositiveIntegerField(default=None, choices=[(1, 'Pinned'), (2, 'Banned')], null=True)),
                ('participate', models.PositiveIntegerField(default=None, choices=[(1, 'Joined'), (2, 'Left')], null=True)),
            ],
            bases=(etalia.core.mixins.ModelDiffMixin, models.Model),
        ),
        migrations.CreateModel(
            name='ThreadUserHistory',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('difference', models.CharField(default='', max_length=256)),
                ('date', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadUserInvite',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('status', models.IntegerField(default=1, choices=[(1, 'Pending'), (2, 'Accepted'), (3, 'Declined')])),
            ],
        ),
    ]
