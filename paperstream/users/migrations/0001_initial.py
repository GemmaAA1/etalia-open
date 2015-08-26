# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from paperstream import users
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
        ('nlp', '0001_initial'),
        ('auth', '0006_require_contenttypes_0002'),
    ]

    operations = [
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('password', models.CharField(max_length=128, verbose_name='password')),
                ('last_login', models.DateTimeField(null=True, blank=True, verbose_name='last login')),
                ('is_superuser', models.BooleanField(verbose_name='superuser status', default=False, help_text='Designates that this user has all permissions without explicitly assigning them.')),
                ('username', models.CharField(blank=True, max_length=255, default='', db_index=True, verbose_name='username (UNUSED)')),
                ('email', models.EmailField(unique=True, max_length=255, verbose_name='Email', db_index=True)),
                ('first_name', models.CharField(blank=True, max_length=255, default='', validators=[users.validators.validate_first_name], verbose_name='First Name')),
                ('last_name', models.CharField(blank=True, max_length=255, default='', validators=[users.validators.validate_last_name], verbose_name='Last Name')),
                ('is_staff', models.BooleanField(verbose_name='staff status', default=False, help_text='Designates whether the user can log into this admin site.')),
                ('is_active', models.BooleanField(verbose_name='active', default=True, help_text='Designates whether this user should be treated as active. Unselect this instead of deleting accounts.')),
            ],
            options={
                'ordering': ('email',),
            },
        ),
        migrations.CreateModel(
            name='Affiliation',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('department', models.CharField(blank=True, max_length=200, default='')),
                ('institution', models.CharField(blank=True, max_length=200, default='')),
                ('city', models.CharField(blank=True, max_length=50, default='')),
                ('state', models.CharField(blank=True, max_length=20, default='')),
                ('country', models.CharField(blank=True, max_length=50, default='')),
            ],
        ),
        migrations.CreateModel(
            name='UserLibJournal',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('papers_in_journal', models.IntegerField(default=0)),
                ('journal', models.ForeignKey(to='library.Journal')),
            ],
        ),
        migrations.CreateModel(
            name='UserLibPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('date_created', models.DateField(null=True, default=None)),
                ('date_last_modified', models.DateField(null=True, default=None)),
                ('authored', models.NullBooleanField(default=None)),
                ('starred', models.NullBooleanField(default=None)),
                ('scored', models.FloatField(default=0.0)),
                ('paper', models.ForeignKey(to='library.Paper', related_name='userlib_paper')),
            ],
            options={
                'ordering': ['-date_created'],
            },
        ),
        migrations.CreateModel(
            name='UserStats',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('state', models.CharField(choices=[('LIN', 'Log in'), ('LOU', 'Log out'), ('LSS', 'Library starts syncing'), ('LES', 'Library ends syncing'), ('FSS', 'Feed starts sync'), ('FES', 'Feed ends sync'), ('EMA', 'Email validated'), ('CRE', 'Create user')], max_length=3)),
                ('options', models.CharField(blank=True, max_length=64, default='')),
                ('datetime', models.DateTimeField(auto_now_add=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserLib',
            fields=[
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL, related_name='lib', serialize=False, primary_key=True)),
                ('state', models.CharField(choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3, default='NON', blank=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='UserSettings',
            fields=[
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL, related_name='settings', serialize=False, primary_key=True)),
                ('scoring_method', models.IntegerField(choices=[(1, 'Average'), (2, 'Average Threshold binary'), (3, 'Average weighted journal'), (4, 'Average weighted journal-date')], default=1, verbose_name='Scoring Algo')),
                ('time_lapse', models.IntegerField(choices=[(61, '2 Months'), (-1, 'All')], default=61, verbose_name='In the past for')),
                ('model', models.ForeignKey(to='nlp.Model', verbose_name='NLP Model')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='userstats',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='stats'),
        ),
        migrations.AlterUniqueTogether(
            name='affiliation',
            unique_together=set([('department', 'institution', 'city', 'state', 'country')]),
        ),
        migrations.AddField(
            model_name='user',
            name='affiliation',
            field=models.ForeignKey(to='users.Affiliation', null=True, default=None),
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(to='auth.Group', blank=True, related_name='user_set', related_query_name='user', help_text='The groups this user belongs to. A user will get all permissions granted to each of their groups.', verbose_name='groups'),
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(to='auth.Permission', blank=True, related_name='user_set', related_query_name='user', help_text='Specific permissions for this user.', verbose_name='user permissions'),
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='userlib',
            field=models.ForeignKey(to='users.UserLib', related_name='userlib_paper'),
        ),
        migrations.AddField(
            model_name='userlibjournal',
            name='userlib',
            field=models.ForeignKey(to='users.UserLib'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='journals',
            field=models.ManyToManyField(to='library.Journal', through='users.UserLibJournal'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='papers',
            field=models.ManyToManyField(to='library.Paper', through='users.UserLibPaper'),
        ),
        migrations.AlterUniqueTogether(
            name='userlibpaper',
            unique_together=set([('userlib', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='userlibjournal',
            unique_together=set([('userlib', 'journal')]),
        ),
    ]
