# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_paper_id_isbn'),
        ('feeds', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='userfeedpaper',
            options={'ordering': ['-score']},
        ),
        migrations.RemoveField(
            model_name='userfeed',
            name='papers',
        ),
        migrations.AddField(
            model_name='userfeed',
            name='name',
            field=models.CharField(default='Main', max_length=100),
        ),
        migrations.AddField(
            model_name='userfeed',
            name='paper_in',
            field=models.ManyToManyField(to='library.Paper', related_name='paper_in'),
        ),
        migrations.AddField(
            model_name='userfeed',
            name='paper_out',
            field=models.ManyToManyField(to='library.Paper', through='feeds.UserFeedPaper', related_name='paper_out'),
        ),
    ]
