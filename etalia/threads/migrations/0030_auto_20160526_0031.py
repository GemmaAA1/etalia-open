# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0010_auto_20160526_0031'),
        ('threads', '0029_auto_20160513_0042'),
    ]

    operations = [
        migrations.CreateModel(
            name='ThreadNeighbors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('time_lapse', models.IntegerField(verbose_name='Days from right now', default=-1, choices=[(7, 'Week'), (30, 'Month'), (60, 'Two Months'), (180, 'Six Months'), (365, 'Year')])),
                ('neighbors', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), blank=True, null=True, size=10)),
                ('ms', models.ForeignKey(to='nlp.MostSimilarThread')),
                ('thread', models.ForeignKey(to='threads.Thread')),
            ],
        ),
        migrations.CreateModel(
            name='ThreadVectors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), null=True, size=300, db_index=True)),
                ('model', models.ForeignKey(to='nlp.Model')),
                ('thread', models.ForeignKey(to='threads.Thread', related_name='vectors')),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='threadvectors',
            unique_together=set([('thread', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadneighbors',
            unique_together=set([('time_lapse', 'thread', 'ms')]),
        ),
    ]
