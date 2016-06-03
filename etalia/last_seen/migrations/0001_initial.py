# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('sites', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='LastSeen',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('module', models.CharField(default='default', max_length=20)),
                ('last_seen', models.DateTimeField(default=django.utils.timezone.now)),
                ('site', models.ForeignKey(to='sites.Site')),
            ],
            options={
                'ordering': ('-last_seen',),
            },
        ),
    ]
