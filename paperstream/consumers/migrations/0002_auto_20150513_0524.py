# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import model_utils.fields
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0016_auto_20150513_0501'),
        ('consumers', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Consumer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True, serialize=False)),
                ('name', models.CharField(max_length=200)),
                ('max_ret', models.IntegerField(default=25)),
                ('day0', models.IntegerField(default=60)),
            ],
        ),
        migrations.CreateModel(
            name='ConsumerJournal',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True, serialize=False)),
                ('status', model_utils.fields.StatusField(default='inactive', no_check_for_status=True, choices=[('inactive', 'inactive'), ('waiting', 'waiting'), ('in_queue', 'in_queue')], max_length=100)),
                ('status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='status')),
            ],
        ),
        migrations.RemoveField(
            model_name='pubmedconsumer',
            name='day0',
        ),
        migrations.RemoveField(
            model_name='pubmedconsumer',
            name='id',
        ),
        migrations.RemoveField(
            model_name='pubmedconsumer',
            name='max_ret',
        ),
        migrations.RemoveField(
            model_name='pubmedconsumer',
            name='name',
        ),
        migrations.CreateModel(
            name='ElsevierConsumer',
            fields=[
                ('consumer_ptr', models.OneToOneField(auto_created=True, serialize=False, parent_link=True, primary_key=True, to='consumers.Consumer')),
                ('api_key', models.CharField(max_length=200)),
            ],
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
        migrations.AddField(
            model_name='pubmedconsumer',
            name='consumer_ptr',
            field=models.OneToOneField(default=1, auto_created=True, serialize=False, parent_link=True, primary_key=True, to='consumers.Consumer'),
            preserve_default=False,
        ),
    ]
