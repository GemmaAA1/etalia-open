# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0014_auto_20150716_1816'),
        ('feeds', '0007_auto_20150713_1829'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserFeedVector',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300)),
                ('feed', models.ForeignKey(to='feeds.UserFeed')),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
