# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0007_auto_20150926_1516'),
    ]

    operations = [
        migrations.CreateModel(
            name='MostSimilarStatus',
            fields=[
                ('id', models.AutoField(serialize=False, verbose_name='ID', primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('host', models.CharField(max_length=512)),
                ('is_downloading', models.BooleanField(default=False)),
                ('ms', models.ForeignKey(to='nlp.MostSimilar')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='mostsimilarstatus',
            unique_together=set([('host', 'ms')]),
        ),
    ]
