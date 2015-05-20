# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Consumer',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('type', models.CharField(max_length=4, choices=[('PUBM', 'PubMed'), ('ELSE', 'Elsevier'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE')])),
                ('name', models.CharField(max_length=200, unique=True)),
                ('ret_max', models.IntegerField(default=25)),
                ('day0', models.IntegerField(default=61)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ConsumerJournal',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('last_date_cons', models.DateTimeField(null=True, default=datetime.datetime(2015, 3, 21, 5, 58, 21, 800053, tzinfo=utc))),
                ('last_number_papers_retrieved', models.IntegerField(default=0)),
                ('last_number_papers_fetched', models.IntegerField(default=0)),
                ('base_countdown_day', models.IntegerField(default=1)),
                ('countdown_day', models.IntegerField(default=0)),
                ('status', model_utils.fields.StatusField(no_check_for_status=True, default='inactive', max_length=100, choices=[('inactive', 'inactive'), ('idle', 'idle'), ('in_queue', 'in_queue'), ('consuming', 'consuming'), ('error', 'error')])),
                ('status_changed', model_utils.fields.MonitorField(monitor='status', default=django.utils.timezone.now)),
            ],
        ),
        migrations.CreateModel(
            name='ConsumerJournalStat',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('date', models.DateTimeField()),
                ('number_papers_fetched', models.IntegerField(default=0)),
                ('number_papers_recorded', models.IntegerField(default=0)),
                ('status', models.CharField(max_length=3, choices=[('SUC', 'Success'), ('FAI', 'Failed')])),
                ('consumer_journal', models.ForeignKey(to='consumers.ConsumerJournal', related_name='stats')),
            ],
        ),
        migrations.CreateModel(
            name='ConsumerArxiv',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, primary_key=True, to='consumers.Consumer', parent_link=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerElsevier',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, primary_key=True, to='consumers.Consumer', parent_link=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerPubmed',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, primary_key=True, to='consumers.Consumer', parent_link=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.AddField(
            model_name='consumerjournal',
            name='consumer',
            field=models.ForeignKey(to='consumers.Consumer'),
        ),
        migrations.AddField(
            model_name='consumerjournal',
            name='journal',
            field=models.ForeignKey(to='library.Journal'),
        ),
        migrations.AddField(
            model_name='consumer',
            name='journals',
            field=models.ManyToManyField(to='library.Journal', through='consumers.ConsumerJournal'),
        ),
        migrations.AlterUniqueTogether(
            name='consumerjournal',
            unique_together=set([('journal', 'consumer')]),
        ),
    ]
