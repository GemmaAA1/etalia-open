# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0023_auto_20160404_2013'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadScore',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0)),
                ('new', models.BooleanField(default=True)),
                ('thread', models.ForeignKey(to='threads.Thread')),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
        migrations.AlterUniqueTogether(
            name='threadfeedthread',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='threadfeedthread',
            name='thread',
        ),
        migrations.RemoveField(
            model_name='threadfeedthread',
            name='thread_feed',
        ),
        migrations.RemoveField(
            model_name='threadfeed',
            name='threads',
        ),
        migrations.DeleteModel(
            name='ThreadFeedThread',
        ),
        migrations.AddField(
            model_name='threadscore',
            name='thread_feed',
            field=models.ForeignKey(to='threads.ThreadFeed'),
        ),
        migrations.AddField(
            model_name='threadfeed',
            name='thread_scores',
            field=models.ManyToManyField(through='threads.ThreadScore', to='threads.Thread'),
        ),
        migrations.AlterUniqueTogether(
            name='threadscore',
            unique_together=set([('thread', 'thread_feed')]),
        ),
    ]
