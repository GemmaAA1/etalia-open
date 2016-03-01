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
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('type', models.CharField(choices=[('PUB', 'PubMed'), ('ELS', 'Elsevier'), ('ARX', 'Arxiv'), ('IEE', 'IEEE')], max_length=3)),
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
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('last_date_cons', models.DateTimeField(default=None, blank=True, null=True)),
                ('last_number_papers_recorded', models.IntegerField(default=0)),
                ('last_number_papers_fetched', models.IntegerField(default=0)),
                ('base_coundown_period', models.IntegerField(default=1)),
                ('coundown_period', models.IntegerField(default=0)),
                ('status', model_utils.fields.StatusField(default='inactive', no_check_for_status=True, choices=[('inactive', 'inactive'), ('idle', 'idle'), ('in_queue', 'in_queue'), ('consuming', 'consuming'), ('error', 'error')], max_length=100)),
                ('status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='status')),
            ],
            options={
                'ordering': ['last_date_cons'],
            },
        ),
        migrations.CreateModel(
            name='ConsumerJournalStat',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('datetime', models.DateTimeField(auto_now_add=True)),
                ('status', models.CharField(choices=[('SUC', 'Success'), ('FAI', 'Failed'), ('RES', 'Reset'), ('ACT', 'Activate'), ('DEA', 'Deactivate')], max_length=3)),
                ('number_papers_fetched', models.IntegerField(default=0)),
                ('number_papers_recorded', models.IntegerField(default=0)),
                ('message', models.CharField(default='', max_length=512)),
                ('consumer_journal', models.ForeignKey(related_name='stats', to='consumers.ConsumerJournal')),
            ],
            options={
                'ordering': ['datetime'],
            },
        ),
        migrations.CreateModel(
            name='ConsumerArxiv',
            fields=[
                ('consumer_ptr', models.OneToOneField(primary_key=True, to='consumers.Consumer', serialize=False, parent_link=True, auto_created=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerElsevier',
            fields=[
                ('consumer_ptr', models.OneToOneField(primary_key=True, to='consumers.Consumer', serialize=False, parent_link=True, auto_created=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.CreateModel(
            name='ConsumerPubmed',
            fields=[
                ('consumer_ptr', models.OneToOneField(primary_key=True, to='consumers.Consumer', serialize=False, parent_link=True, auto_created=True)),
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
        migrations.AlterUniqueTogether(
            name='consumerjournal',
            unique_together=set([('journal', 'consumer')]),
        ),
    ]
