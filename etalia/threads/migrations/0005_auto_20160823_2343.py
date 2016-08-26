# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0004_auto_20160605_2237'),
    ]

    operations = [
        migrations.CreateModel(
            name='PubPeer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('doi', models.CharField(db_index=True, default='', blank=True, max_length=64, verbose_name='DOI', unique=True)),
                ('link', models.URLField()),
                ('pubpeer_id', models.CharField(max_length=30)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='PubPeerComments',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('body', models.TextField()),
                ('date', models.DateTimeField()),
                ('pubpeercomment_id', models.IntegerField()),
                ('permalink', models.URLField()),
                ('rating', models.FloatField()),
                ('user', models.CharField(max_length=128)),
                ('pubpeer', models.ForeignKey(related_name='comments', to='threads.PubPeer')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='thread',
            name='pubpeer',
            field=models.ForeignKey(to='threads.PubPeer', blank=True, null=True),
        ),
    ]
