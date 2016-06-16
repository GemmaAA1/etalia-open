# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20160605_1941'),
    ]

    operations = [
        migrations.CreateModel(
            name='PaperUserHistory',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('difference', models.CharField(max_length=256, default='')),
                ('date', models.DateTimeField(auto_now_add=True)),
                ('paperuser', models.ForeignKey(to='library.PaperUser', related_name='history')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
