# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0026_auto_20160421_0046'),
        ('nlp', '0005_auto_20151208_0359'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadNeighbor',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')], verbose_name='Days from right now')),
                ('neighbors', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), size=10, null=True, blank=True)),
                ('thread', models.ForeignKey(to='threads.Thread')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ThreadVectors',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), db_index=True, size=300, null=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('thread', models.ForeignKey(related_name='vectors', to='threads.Thread')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='threadvectors',
            unique_together=set([('thread', 'model')]),
        ),
    ]
