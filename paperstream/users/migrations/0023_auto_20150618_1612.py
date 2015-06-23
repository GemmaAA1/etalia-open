# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0022_auto_20150616_2232'),
    ]

    operations = [
        migrations.AddField(
            model_name='affiliation',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 19, 178032, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='affiliation',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 23, 393815, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userlib',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 25, 536264, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userlib',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 27, 224267, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userlibjournal',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 28, 800239, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userlibjournal',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 32, 448211, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 34, 17527, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='userlibpaper',
            name='modified',
            field=models.DateTimeField(default=datetime.datetime(2015, 6, 18, 16, 12, 36, 160463, tzinfo=utc), auto_now=True),
            preserve_default=False,
        ),
    ]
