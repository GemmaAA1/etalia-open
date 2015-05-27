# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
        ('users', '0007_auto_20150522_1937'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserLib',
            fields=[
                ('user', models.OneToOneField(related_name='profile', serialize=False, to=settings.AUTH_USER_MODEL, primary_key=True)),
                ('library_status', models.CharField(default='', choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')], blank=True, max_length=3)),
                ('feed_status', models.CharField(default='', choices=[('', 'Uninitialized'), ('SYN', 'Synced'), ('ING', 'Syncing')], blank=True, max_length=3)),
                ('is_library_hooked', models.BooleanField(default=False)),
                ('affiliation', models.ForeignKey(to='users.Affiliation', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserLibJournal',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('papers_in_journal', models.IntegerField(default=0)),
                ('score', models.FloatField(default=0.0)),
                ('journal', models.ForeignKey(to='library.Journal')),
                ('user', models.ForeignKey(to='users.UserLib')),
            ],
        ),
        migrations.CreateModel(
            name='UserLibPaper',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('date_added', models.DateField(default=datetime.date(2000, 1, 1))),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('user', models.ForeignKey(to='users.UserLib')),
            ],
            options={
                'ordering': ['-date_added'],
            },
        ),
        migrations.CreateModel(
            name='UserLibStats',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', auto_created=True, primary_key=True)),
                ('state', models.CharField(choices=[('LOG', 'Log in'), ('LIB', 'Library syncing'), ('FEE', 'Feed syncing')], max_length=3)),
                ('datetime', models.DateTimeField(auto_now_add=True)),
                ('user', models.ForeignKey(related_name='stats', to='users.UserLib')),
            ],
        ),
        migrations.RemoveField(
            model_name='userjournal',
            name='journal',
        ),
        migrations.RemoveField(
            model_name='userjournal',
            name='user',
        ),
        migrations.RemoveField(
            model_name='userpaper',
            name='paper',
        ),
        migrations.RemoveField(
            model_name='userpaper',
            name='user',
        ),
        migrations.RemoveField(
            model_name='userstats',
            name='user',
        ),
        migrations.RemoveField(
            model_name='user',
            name='affiliation',
        ),
        migrations.RemoveField(
            model_name='user',
            name='feed_status',
        ),
        migrations.RemoveField(
            model_name='user',
            name='is_library_hooked',
        ),
        migrations.RemoveField(
            model_name='user',
            name='journals',
        ),
        migrations.RemoveField(
            model_name='user',
            name='library_status',
        ),
        migrations.RemoveField(
            model_name='user',
            name='own_papers',
        ),
        migrations.RemoveField(
            model_name='user',
            name='papers',
        ),
        migrations.DeleteModel(
            name='UserJournal',
        ),
        migrations.DeleteModel(
            name='UserPaper',
        ),
        migrations.DeleteModel(
            name='UserStats',
        ),
        migrations.AddField(
            model_name='userlib',
            name='journals',
            field=models.ManyToManyField(to='library.Journal', through='users.UserLibJournal'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='own_papers',
            field=models.ManyToManyField(related_name='own', to='library.Paper'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='users.UserLibPaper'),
        ),
    ]
