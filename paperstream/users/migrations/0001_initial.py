# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import jsonfield.fields
import paperstream.users.validators
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('auth', '0006_require_contenttypes_0002'),
        ('nlp', '0001_initial'),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('password', models.CharField(verbose_name='password', max_length=128)),
                ('last_login', models.DateTimeField(blank=True, null=True, verbose_name='last login')),
                ('is_superuser', models.BooleanField(default=False, verbose_name='superuser status', help_text='Designates that this user has all permissions without explicitly assigning them.')),
                ('username', models.CharField(default='', blank=True, verbose_name='username (UNUSED)', max_length=255, db_index=True)),
                ('email', models.EmailField(unique=True, verbose_name='Email', max_length=255, db_index=True)),
                ('first_name', models.CharField(default='', blank=True, verbose_name='First Name', validators=[paperstream.users.validators.validate_first_name], max_length=255)),
                ('last_name', models.CharField(default='', blank=True, verbose_name='Last Name', validators=[paperstream.users.validators.validate_last_name], max_length=255)),
                ('title', models.CharField(default='', blank=True, max_length=32)),
                ('position', models.CharField(default='', blank=True, max_length=64)),
                ('is_staff', models.BooleanField(default=False, verbose_name='staff status', help_text='Designates whether the user can log into this admin site.')),
                ('is_active', models.BooleanField(default=True, verbose_name='active', help_text='Designates whether this user should be treated as active. Unselect this instead of deleting accounts.')),
                ('init_step', models.CharField(default='NON', choices=[('NON', 'uninitialized'), ('LIB', 'library'), ('STR', 'personalized stream'), ('TRE', 'trends'), ('IDL', 'done')], max_length=3, help_text='Tag where init user stands')),
                ('photo', models.ImageField(null=True, upload_to='photos')),
            ],
            options={
                'ordering': ('email',),
            },
        ),
        migrations.CreateModel(
            name='Affiliation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('department', models.CharField(default='', blank=True, max_length=200)),
                ('institution', models.CharField(default='', blank=True, max_length=200)),
                ('city', models.CharField(default='', blank=True, max_length=50)),
                ('state', models.CharField(default='', blank=True, max_length=20)),
                ('country', models.CharField(default='', blank=True, max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='FeedLayout',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('stream_filter', jsonfield.fields.JSONField(null=True)),
                ('trend_filter', jsonfield.fields.JSONField(null=True)),
                ('library_filter', jsonfield.fields.JSONField(null=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='UserLibJournal',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('papers_in_journal', models.IntegerField(default=0)),
                ('journal', models.ForeignKey(to='library.Journal')),
            ],
            options={
                'ordering': ('-papers_in_journal',),
            },
        ),
        migrations.CreateModel(
            name='UserLibPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('date_created', models.DateField(default=None, null=True)),
                ('date_last_modified', models.DateField(default=None, null=True)),
                ('authored', models.NullBooleanField(default=None)),
                ('starred', models.NullBooleanField(default=None)),
                ('scored', models.FloatField(default=0.0)),
                ('paper_provider_id', models.CharField(default='', max_length=64)),
                ('is_trashed', models.BooleanField(default=False)),
                ('paper', models.ForeignKey(related_name='userlib_paper', to='library.Paper')),
            ],
            options={
                'ordering': ['-date_created'],
            },
        ),
        migrations.CreateModel(
            name='UserStats',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('message', models.CharField(max_length=128)),
                ('options', models.CharField(default='', blank=True, max_length=128)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='UserTaste',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('context_source', models.CharField(max_length=128)),
                ('is_ticked', models.BooleanField(default=False)),
                ('is_liked', models.BooleanField(default=False)),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='UserLib',
            fields=[
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('user', models.OneToOneField(primary_key=True, related_name='lib', serialize=False, to=settings.AUTH_USER_MODEL)),
                ('state', models.CharField(default='NON', blank=True, choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3)),
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
                ('user', models.OneToOneField(primary_key=True, related_name='settings', serialize=False, to=settings.AUTH_USER_MODEL)),
                ('stream_method', models.IntegerField(default=1, verbose_name='Method')),
                ('stream_time_lapse', models.IntegerField(default=30, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (365, 'Year'), (-1, 'All')], verbose_name='In the past for')),
                ('trend_method', models.IntegerField(default=1, verbose_name='Method')),
                ('trend_time_lapse', models.IntegerField(default=30, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (365, 'Year'), (-1, 'All')], verbose_name='In the past for')),
                ('stream_model', models.ForeignKey(related_name='stream_model', verbose_name='NLP Model', to='nlp.Model')),
                ('trend_model', models.ForeignKey(related_name='trend_model', verbose_name='NLP Model', to='nlp.Model')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='usertaste',
            name='user',
            field=models.ForeignKey(related_name='tastes', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='userstats',
            name='user',
            field=models.ForeignKey(related_name='stats', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='feedlayout',
            name='user',
            field=models.OneToOneField(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AlterUniqueTogether(
            name='affiliation',
            unique_together=set([('department', 'institution', 'city', 'state', 'country')]),
        ),
        migrations.AddField(
            model_name='user',
            name='affiliation',
            field=models.ForeignKey(default=None, to='users.Affiliation', null=True),
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(blank=True, related_name='user_set', related_query_name='user', to='auth.Group', verbose_name='groups', help_text='The groups this user belongs to. A user will get all permissions granted to each of their groups.'),
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(blank=True, related_name='user_set', related_query_name='user', to='auth.Permission', verbose_name='user permissions', help_text='Specific permissions for this user.'),
        ),
        migrations.AlterUniqueTogether(
            name='usertaste',
            unique_together=set([('user', 'paper')]),
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='userlib',
            field=models.ForeignKey(related_name='userlib_paper', to='users.UserLib'),
        ),
        migrations.AddField(
            model_name='userlibjournal',
            name='userlib',
            field=models.ForeignKey(to='users.UserLib'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='journals',
            field=models.ManyToManyField(through='users.UserLibJournal', to='library.Journal'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='papers',
            field=models.ManyToManyField(through='users.UserLibPaper', to='library.Paper'),
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
