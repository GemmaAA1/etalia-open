# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_auto_20150929_0710'),
    ]

    operations = [
        migrations.CreateModel(
            name='Stats',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True, serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('nb_papers', models.IntegerField(default=0)),
                ('nb_journal', models.IntegerField(default=0)),
                ('nb_authors', models.IntegerField(default=0)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
