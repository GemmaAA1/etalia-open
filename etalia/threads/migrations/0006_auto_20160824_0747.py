# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0005_auto_20160823_2343'),
    ]

    operations = [
        migrations.CreateModel(
            name='PubPeerComment',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('body', models.TextField(default='')),
                ('date', models.DateTimeField()),
                ('pubpeercomment_id', models.IntegerField(unique=True)),
                ('permalink', models.URLField()),
                ('rating', models.FloatField()),
                ('user', models.CharField(max_length=128)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='pubpeercomments',
            name='pubpeer',
        ),
        migrations.RemoveField(
            model_name='thread',
            name='pubpeer',
        ),
        migrations.AddField(
            model_name='pubpeer',
            name='thread',
            field=models.OneToOneField(default=1, to='threads.Thread'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='pubpeer',
            name='pubpeer_id',
            field=models.CharField(max_length=30, unique=True),
        ),
        migrations.DeleteModel(
            name='PubPeerComments',
        ),
        migrations.AddField(
            model_name='pubpeercomment',
            name='pubpeer',
            field=models.ForeignKey(related_name='comments', to='threads.PubPeer'),
        ),
    ]
