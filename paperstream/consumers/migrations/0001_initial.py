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
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('type', models.CharField(max_length=3, choices=[('PUB', 'PubMed'), ('ELS', 'Elsevier'), ('ARX', 'Arxiv'), ('IEE', 'IEEE')])),
                ('name', models.CharField(max_length=200, unique=True)),
                ('ret_max', models.IntegerField(default=25)),
                ('day0', models.IntegerField(default=365)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ConsumerJournal',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('last_date_cons', models.DateTimeField(null=True, blank=True)),
                ('last_number_papers_recorded', models.IntegerField(default=0)),
                ('last_number_papers_fetched', models.IntegerField(default=0)),
                ('base_coundown_period', models.IntegerField(default=1)),
                ('coundown_period', models.IntegerField(default=0)),
                ('status', model_utils.fields.StatusField(default='inactive', no_check_for_status=True, max_length=100, choices=[('inactive', 'inactive'), ('idle', 'idle'), ('in_queue', 'in_queue'), ('consuming', 'consuming'), ('error', 'error')])),
                ('status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='status')),
            ],
            options={
                'ordering': ['last_date_cons'],
            },
        ),
        migrations.CreateModel(
            name='ConsumerJournalStat',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, serialize=False, verbose_name='ID')),
                ('datetime', models.DateTimeField(auto_now_add=True)),
                ('number_papers_fetched', models.IntegerField()),
                ('number_papers_recorded', models.IntegerField()),
                ('status', models.CharField(max_length=3, choices=[('SUC', 'Success'), ('FAI', 'Failed')])),
                ('consumer_journal', models.ForeignKey(related_name='stats', to='consumers.ConsumerJournal')),
            ],
        ),
        migrations.CreateModel(
            name='ConsumerArxiv',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, to='consumers.Consumer', primary_key=True, parent_link=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerElsevier',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, to='consumers.Consumer', primary_key=True, parent_link=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerPubmed',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, to='consumers.Consumer', primary_key=True, parent_link=True)),
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
            field=models.ManyToManyField(through='consumers.ConsumerJournal', to='library.Journal'),
        ),
    ]
