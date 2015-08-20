# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0006_auto_20150710_0053'),
        ('feeds', '0015_auto_20150817_1717'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserFeedMatchPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('is_score_computed', models.BooleanField(default=False)),
                ('is_disliked', models.BooleanField(default=False)),
                ('is_liked', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.AlterUniqueTogether(
            name='userfeedpaper',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='userfeedpaper',
            name='feed',
        ),
        migrations.RemoveField(
            model_name='userfeedpaper',
            name='paper',
        ),
        migrations.AlterField(
            model_name='userfeed',
            name='papers_match',
            field=models.ManyToManyField(to='library.Paper', related_name='feed_match', through='feeds.UserFeedMatchPaper'),
        ),
        migrations.DeleteModel(
            name='UserFeedPaper',
        ),
        migrations.AddField(
            model_name='userfeedmatchpaper',
            name='feed',
            field=models.ForeignKey(to='feeds.UserFeed'),
        ),
        migrations.AddField(
            model_name='userfeedmatchpaper',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AlterUniqueTogether(
            name='userfeedmatchpaper',
            unique_together=set([('feed', 'paper')]),
        ),
    ]
