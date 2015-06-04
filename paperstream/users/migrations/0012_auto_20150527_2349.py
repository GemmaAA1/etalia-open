# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0011_user_username'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='userlibpaper',
            options={'ordering': ['-date_created']},
        ),
        migrations.RenameField(
            model_name='userlibpaper',
            old_name='date_added',
            new_name='date_created',
        ),
        migrations.RemoveField(
            model_name='userlib',
            name='own_papers',
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='authored',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='date_last_modified',
            field=models.DateField(default=datetime.date(2000, 1, 1)),
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='scored',
            field=models.FloatField(default=0.0),
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='starred',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='userlibstats',
            name='number_papers',
            field=models.IntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='userlibstats',
            name='state',
            field=models.CharField(max_length=3, choices=[('LOG', 'Log in'), ('LIB', 'Library sync'), ('FEE', 'Feed sync')]),
        ),
        migrations.AlterUniqueTogether(
            name='affiliation',
            unique_together=set([('department', 'institution', 'city', 'state', 'country')]),
        ),
    ]
