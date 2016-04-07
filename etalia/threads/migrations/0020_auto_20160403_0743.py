# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0019_auto_20160402_2354'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadUserHistory',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('date', models.DateTimeField(auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('action', models.CharField(max_length=128)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AlterModelOptions(
            name='threaduser',
            options={},
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='first_joined_at',
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='is_banned',
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='is_joined',
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='is_left',
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='is_pinned',
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='last_left_at',
        ),
        migrations.RemoveField(
            model_name='threaduser',
            name='num_comments',
        ),
        migrations.AddField(
            model_name='threaduser',
            name='participate',
            field=models.PositiveIntegerField(null=True, default=None, choices=[(1, 'Joined'), (2, 'Left')]),
        ),
        migrations.AddField(
            model_name='threaduser',
            name='watch',
            field=models.PositiveIntegerField(null=True, default=None, choices=[(1, 'Pinned'), (2, 'Banned')]),
        ),
        migrations.AddField(
            model_name='threaduserhistory',
            name='threaduser',
            field=models.ForeignKey(to='threads.ThreadUser'),
        ),
    ]
