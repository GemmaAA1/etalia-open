# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0018_auto_20150606_0412'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userstats',
            name='number_papers',
            field=models.IntegerField(blank=True, default=None, null=True),
        ),
        migrations.AlterField(
            model_name='userstats',
            name='state',
            field=models.CharField(choices=[('LIN', 'Log in'), ('LOU', 'Log out'), ('LIB', 'Library sync'), ('FEE', 'Feed sync'), ('EMA', 'Email validated'), ('CRE', 'Create user')], max_length=3),
        ),
    ]
