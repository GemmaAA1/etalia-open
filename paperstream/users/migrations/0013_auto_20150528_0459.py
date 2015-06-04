# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0012_auto_20150527_2349'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserStats',
            fields=[
                ('user', models.OneToOneField(serialize=False, to=settings.AUTH_USER_MODEL, primary_key=True, related_name='stats')),
                ('state', models.CharField(max_length=3, choices=[('LIN', 'Log in'), ('LOU', 'Log out'), ('LIB', 'Library sync'), ('FEE', 'Feed sync')])),
                ('number_papers', models.IntegerField(default=0)),
                ('datetime', models.DateTimeField(auto_now_add=True)),
            ],
        ),
        migrations.RemoveField(
            model_name='userlibstats',
            name='userlib',
        ),
        migrations.RemoveField(
            model_name='userlib',
            name='affiliation',
        ),
        migrations.RemoveField(
            model_name='userlib',
            name='is_library_hooked',
        ),
        migrations.AddField(
            model_name='user',
            name='affiliation',
            field=models.ForeignKey(null=True, to='users.Affiliation', default=None),
        ),
        migrations.AlterField(
            model_name='userlib',
            name='feed_status',
            field=models.CharField(max_length=3, choices=[('', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], blank=True, default=''),
        ),
        migrations.AlterField(
            model_name='userlib',
            name='library_status',
            field=models.CharField(max_length=3, choices=[('', 'Uninitialized'), ('IDL', 'Idle'), ('ING', 'Syncing')], blank=True, default=''),
        ),
        migrations.AlterField(
            model_name='userlib',
            name='user',
            field=models.OneToOneField(serialize=False, to=settings.AUTH_USER_MODEL, primary_key=True, related_name='lib'),
        ),
        migrations.DeleteModel(
            name='UserLibStats',
        ),
    ]
