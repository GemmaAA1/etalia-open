# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('auth', '0001_initial'),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('password', models.CharField(verbose_name='password', max_length=128)),
                ('last_login', models.DateTimeField(verbose_name='last login', null=True, blank=True)),
                ('is_superuser', models.BooleanField(verbose_name='superuser status', default=False, help_text='Designates that this user has all permissions without explicitly assigning them.')),
                ('first_name', models.CharField(default='', max_length=200, blank=True)),
                ('last_name', models.CharField(default='', max_length=200, blank=True)),
                ('email', models.EmailField(verbose_name='Email address', max_length=255, unique=True)),
                ('lib_size', models.IntegerField(default=0)),
                ('is_active', models.BooleanField(default=True)),
                ('is_admin', models.BooleanField(default=False)),
                ('library_status', models.CharField(choices=[('SYN', 'Synced'), ('ING', 'Syncing')], max_length=3)),
                ('feed_status', models.CharField(choices=[('SYN', 'Synced'), ('ING', 'Syncing')], max_length=3)),
                ('is_library_hooked', models.BooleanField(default=False)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Affiliation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('department', models.CharField(max_length=200, null=True)),
                ('institution', models.CharField(max_length=200, null=True)),
                ('city', models.CharField(max_length=50, null=True)),
                ('state', models.CharField(max_length=20, null=True)),
                ('country', models.CharField(max_length=50, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserJournal',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('count', models.FloatField(default=0.0)),
                ('last_addition', models.DateField(default=datetime.date(1900, 1, 1))),
                ('score', models.FloatField(default=1.0)),
                ('journal', models.ForeignKey(to='library.Journal')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='UserPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('date_added', models.DateField(default=datetime.date(2000, 1, 1))),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-date_added'],
            },
        ),
        migrations.CreateModel(
            name='UserStats',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('doing', models.CharField(choices=[('LOG', 'Log in'), ('LIB', 'Library syncing'), ('FEE', 'Feed syncing')], max_length=3)),
                ('when', models.DateTimeField()),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='stats')),
            ],
        ),
        migrations.AddField(
            model_name='user',
            name='affiliation',
            field=models.ForeignKey(null=True, to='users.Affiliation'),
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(related_name='user_set', blank=True, verbose_name='groups', help_text='The groups this user belongs to. A user will get all permissions granted to each of their groups.', to='auth.Group', related_query_name='user'),
        ),
        migrations.AddField(
            model_name='user',
            name='journals',
            field=models.ManyToManyField(to='library.Journal', through='users.UserJournal'),
        ),
        migrations.AddField(
            model_name='user',
            name='own_papers',
            field=models.ManyToManyField(related_name='own', to='library.Paper'),
        ),
        migrations.AddField(
            model_name='user',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='users.UserPaper'),
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(related_name='user_set', blank=True, verbose_name='user permissions', help_text='Specific permissions for this user.', to='auth.Permission', related_query_name='user'),
        ),
    ]
