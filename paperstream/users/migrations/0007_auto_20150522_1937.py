# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0006_auto_20150521_2314'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='user',
            options={'ordering': ('email',)},
        ),
        migrations.RemoveField(
            model_name='user',
            name='is_admin',
        ),
        migrations.RemoveField(
            model_name='user',
            name='lib_size',
        ),
        migrations.RemoveField(
            model_name='userjournal',
            name='count',
        ),
        migrations.RemoveField(
            model_name='userjournal',
            name='last_addition',
        ),
        migrations.AddField(
            model_name='user',
            name='is_staff',
            field=models.BooleanField(help_text='Designates whether the user can log into this admin site.', verbose_name='staff status', default=False),
        ),
        migrations.AddField(
            model_name='userjournal',
            name='papers_in_journal',
            field=models.IntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='user',
            name='email',
            field=models.EmailField(max_length=255, db_index=True, verbose_name='email address', unique=True),
        ),
        migrations.AlterField(
            model_name='user',
            name='is_active',
            field=models.BooleanField(help_text='Designates whether this user should be treated as active. Unselect this instead of deleting accounts.', verbose_name='active', default=True),
        ),
        migrations.AlterField(
            model_name='userjournal',
            name='score',
            field=models.FloatField(default=0.0),
        ),
    ]
