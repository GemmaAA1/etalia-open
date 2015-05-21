# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Consumer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('type', models.CharField(choices=[('PUB', 'PubMed'), ('ELS', 'Elsevier'), ('ARX', 'Arxiv'), ('IEE', 'IEEE')], max_length=3)),
                ('name', models.CharField(max_length=200, unique=True)),
                ('ret_max', models.IntegerField(default=25)),
                ('day0', models.IntegerField(default=60)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ConsumerJournal',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('last_date_cons', models.DateTimeField(null=True, blank=True)),
                ('last_number_papers_recorded', models.IntegerField(default=0)),
                ('last_number_papers_fetched', models.IntegerField(default=0)),
                ('base_countdown_day', models.IntegerField(default=1)),
                ('countdown_day', models.IntegerField(default=0)),
                ('status', model_utils.fields.StatusField(choices=[('inactive', 'inactive'), ('idle', 'idle'), ('in_queue', 'in_queue'), ('consuming', 'consuming'), ('error', 'error')], default='inactive', max_length=100, no_check_for_status=True)),
                ('status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='status')),
            ],
        ),
        migrations.CreateModel(
            name='ConsumerJournalStat',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('date', models.DateTimeField()),
                ('number_papers_fetched', models.IntegerField(default=0)),
                ('number_papers_recorded', models.IntegerField(default=0)),
                ('status', models.CharField(choices=[('SUC', 'Success'), ('FAI', 'Failed')], max_length=3)),
                ('consumer_journal', models.ForeignKey(to='consumers.ConsumerJournal', related_name='stats')),
            ],
        ),
        migrations.CreateModel(
            name='ConsumerArxiv',
            fields=[
                ('consumer_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='consumers.Consumer')),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerElsevier',
            fields=[
                ('consumer_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='consumers.Consumer')),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerPubmed',
            fields=[
                ('consumer_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='consumers.Consumer')),
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
