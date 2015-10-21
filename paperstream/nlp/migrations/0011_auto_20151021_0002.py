# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0007_auto_20151009_1916'),
        ('nlp', '0010_auto_20151013_2112'),
    ]

    operations = [
        migrations.CreateModel(
            name='JournalNeighbors',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('neighbors', django.contrib.postgres.fields.ArrayField(blank=True, null=True, base_field=models.IntegerField(null=True), size=10)),
                ('journal', models.ForeignKey(to='library.Journal', related_name='neighbors')),
                ('model', models.ForeignKey(to='nlp.Model')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([('journal', 'model')]),
        ),
    ]
