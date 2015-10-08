# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_stats'),
        ('users', '0005_auto_20151001_0754'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserTaste',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('is_disliked', models.BooleanField(default=False)),
                ('is_liked', models.BooleanField(default=False)),
                ('paper', models.ForeignKey(to='library.Paper')),
                ('user', models.ForeignKey(related_name='tastes', to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='usertaste',
            unique_together=set([('user', 'paper')]),
        ),
    ]
