# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150710_0053'),
        ('feeds', '0016_auto_20150819_1836'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserFeedSeedPaper',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('feed', models.ForeignKey(to='feeds.UserFeed')),
                ('paper', models.ForeignKey(to='library.Paper')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='userfeedseedpaper',
            unique_together=set([('feed', 'paper')]),
        ),
    ]
