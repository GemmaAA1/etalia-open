# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0003_auto_20150513_0527'),
    ]

    operations = [
        migrations.CreateModel(
            name='ArxivConsumer',
            fields=[
                ('consumer_ptr', models.OneToOneField(to='consumers.Consumer', parent_link=True, serialize=False, primary_key=True, auto_created=True)),
                ('url', models.URLField()),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.RemoveField(
            model_name='elsevierconsumer',
            name='api_key',
        ),
        migrations.RemoveField(
            model_name='pubmedconsumer',
            name='email',
        ),
    ]
