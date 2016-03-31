# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0023_auto_20160308_0835'),
    ]

    operations = [
        migrations.CreateModel(
            name='Relationship',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('status', models.IntegerField(default=1, choices=[(1, 'Following'), (2, 'Blocked')])),
                ('from_user', models.ForeignKey(related_name='relation_from_users', to=settings.AUTH_USER_MODEL)),
                ('to_user', models.ForeignKey(related_name='relation_to_users', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_roll_back_deltatime',
            field=models.IntegerField(default=36, verbose_name='Roll-back time (months)'),
        ),
        migrations.AddField(
            model_name='user',
            name='relationships',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL, through='users.Relationship', related_name='related_to'),
        ),
    ]
