# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings
import etalia.users.validators


class Migration(migrations.Migration):

    dependencies = [
        ('auth', '0006_require_contenttypes_0002'),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('password', models.CharField(verbose_name='password', max_length=128)),
                ('last_login', models.DateTimeField(blank=True, verbose_name='last login', null=True)),
                ('is_superuser', models.BooleanField(default=False, help_text='Designates that this user has all permissions without explicitly assigning them.', verbose_name='superuser status')),
                ('username', models.CharField(blank=True, max_length=255, default='', verbose_name='username (UNUSED)', db_index=True)),
                ('email', models.EmailField(max_length=255, unique=True, verbose_name='Email', db_index=True)),
                ('first_name', models.CharField(blank=True, default='', verbose_name='First Name', validators=[etalia.users.validators.validate_first_name], max_length=255)),
                ('last_name', models.CharField(blank=True, default='', verbose_name='Last Name', validators=[etalia.users.validators.validate_last_name], max_length=255)),
                ('title', models.CharField(blank=True, default='', max_length=32)),
                ('position', models.CharField(blank=True, default='', max_length=64)),
                ('is_staff', models.BooleanField(default=False, help_text='Designates whether the user can log into this admin site.', verbose_name='staff status')),
                ('is_active', models.BooleanField(default=True, help_text='Designates whether this user should be treated as active. Unselect this instead of deleting accounts.', verbose_name='active')),
                ('is_alpha', models.BooleanField(default=True, help_text='Designates whether this user should be treated as an early adopter user.', verbose_name='alpha')),
                ('init_step', models.CharField(default='NON', help_text='Tag where init user stands', choices=[('NON', 'uninitialized'), ('LIB', 'library'), ('STR', 'personalized stream'), ('TRE', 'trends'), ('IDL', 'done')], max_length=3)),
            ],
            options={
                'ordering': ('email',),
            },
        ),
        migrations.CreateModel(
            name='Affiliation',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('department', models.CharField(blank=True, default='', max_length=200)),
                ('institution', models.CharField(blank=True, default='', max_length=200)),
                ('city', models.CharField(blank=True, default='', max_length=50)),
                ('state', models.CharField(blank=True, default='', max_length=20)),
                ('country', models.CharField(blank=True, default='', max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Relationship',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('status', models.IntegerField(default=1, choices=[(1, 'Following'), (2, 'Blocked')])),
            ],
        ),
        migrations.CreateModel(
            name='UserLibAuthor',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('occurrence', models.IntegerField(default=0)),
                ('author', models.ForeignKey(to='library.Author')),
            ],
            options={
                'ordering': ('-occurrence',),
            },
        ),
        migrations.CreateModel(
            name='UserLibJournal',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('occurrence', models.IntegerField(default=0)),
                ('journal', models.ForeignKey(to='library.Journal')),
            ],
            options={
                'ordering': ('-occurrence',),
            },
        ),
        migrations.CreateModel(
            name='UserLibPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('date_created', models.DateField(default=None, null=True)),
                ('date_last_modified', models.DateField(default=None, null=True)),
                ('authored', models.NullBooleanField(default=None)),
                ('starred', models.NullBooleanField(default=None)),
                ('scored', models.FloatField(default=0.0)),
                ('paper_provider_id', models.CharField(default='', max_length=64)),
                ('is_trashed', models.BooleanField(default=False)),
                ('paper', models.ForeignKey(to='library.Paper', related_name='userlib_paper')),
            ],
            options={
                'ordering': ['-date_created'],
            },
        ),
        migrations.CreateModel(
            name='UserStats',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('message', models.CharField(max_length=128)),
                ('options', models.CharField(blank=True, default='', max_length=128)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='UserTaste',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('source', models.CharField(max_length=128)),
                ('is_banned', models.BooleanField(default=False)),
                ('is_pinned', models.BooleanField(default=False)),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='UserLib',
            fields=[
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL, primary_key=True, related_name='lib', serialize=False)),
                ('state', models.CharField(blank=True, default='NON', choices=[('NON', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], max_length=3)),
                ('d_oldest', models.DateField(blank=True, null=True)),
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
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL, primary_key=True, related_name='settings', serialize=False)),
                ('stream_method', models.IntegerField(default=0, choices=[(0, 'Content Based Scoring (simple)')], verbose_name='Method')),
                ('stream_author_weight', models.FloatField(default=1.0, verbose_name='Author weight')),
                ('stream_journal_weight', models.FloatField(default=1.0, verbose_name='Journal weight')),
                ('stream_vector_weight', models.FloatField(default=1.0, verbose_name='Content weight')),
                ('stream_roll_back_deltatime', models.IntegerField(default=36, verbose_name='Roll-back time (months)')),
                ('trend_method', models.IntegerField(default=0, choices=[(0, 'Method #1')], verbose_name='Method')),
                ('trend_doc_weight', models.FloatField(default=1.0, verbose_name='Content weight')),
                ('trend_altmetric_weight', models.FloatField(default=1.0, verbose_name='Altmetric weight')),
                ('email_digest_frequency', models.IntegerField(default=7, choices=[(7, 'Weekly'), (15, 'Bi-Weekly'), (30, 'Monthly'), (30, 'Bi-Monthly'), (-1, 'Never')], verbose_name='Email digest frequency')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='usertaste',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='tastes'),
        ),
        migrations.AddField(
            model_name='userstats',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='stats'),
        ),
        migrations.AddField(
            model_name='relationship',
            name='from_user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='relation_from_users'),
        ),
        migrations.AddField(
            model_name='relationship',
            name='to_user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='relation_to_users'),
        ),
        migrations.AlterUniqueTogether(
            name='affiliation',
            unique_together=set([('department', 'institution', 'city', 'state', 'country')]),
        ),
        migrations.AddField(
            model_name='user',
            name='affiliation',
            field=models.ForeignKey(to='users.Affiliation', default=None, null=True),
        ),
        migrations.AddField(
            model_name='user',
            name='groups',
            field=models.ManyToManyField(blank=True, to='auth.Group', help_text='The groups this user belongs to. A user will get all permissions granted to each of their groups.', related_query_name='user', related_name='user_set', verbose_name='groups'),
        ),
        migrations.AddField(
            model_name='user',
            name='relationships',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL, through='users.Relationship', related_name='related_to'),
        ),
        migrations.AddField(
            model_name='user',
            name='user_permissions',
            field=models.ManyToManyField(blank=True, to='auth.Permission', help_text='Specific permissions for this user.', related_query_name='user', related_name='user_set', verbose_name='user permissions'),
        ),
        migrations.AlterUniqueTogether(
            name='usertaste',
            unique_together=set([('user', 'paper')]),
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
            model_name='userlibauthor',
            name='userlib',
            field=models.ForeignKey(to='users.UserLib'),
        ),
        migrations.AddField(
            model_name='userlib',
            name='authors',
            field=models.ManyToManyField(to='library.Author', through='users.UserLibAuthor'),
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
            name='relationship',
            unique_together=set([('from_user', 'to_user')]),
        ),
        migrations.AlterUniqueTogether(
            name='userlibpaper',
            unique_together=set([('userlib', 'paper')]),
        ),
        migrations.AlterUniqueTogether(
            name='userlibjournal',
            unique_together=set([('userlib', 'journal')]),
        ),
        migrations.AlterUniqueTogether(
            name='userlibauthor',
            unique_together=set([('userlib', 'author')]),
        ),
    ]
