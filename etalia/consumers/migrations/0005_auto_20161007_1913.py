# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import etalia.consumers.mixins


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0004_pubpeerconsumer'),
    ]

    operations = [
        migrations.CreateModel(
            name='ConsumerBiorxiv',
            fields=[
                ('consumer_ptr', models.OneToOneField(primary_key=True, serialize=False, to='consumers.Consumer', parent_link=True, auto_created=True)),
            ],
            options={
                'abstract': False,
            },
            bases=(etalia.consumers.mixins.SimpleCrawlerListMixin, 'consumers.consumer'),
        ),
        migrations.AlterField(
            model_name='consumer',
            name='type',
            field=models.CharField(choices=[('PUB', 'PubMed'), ('ELS', 'Elsevier'), ('ARX', 'Arxiv'), ('IEE', 'IEEE'), ('BIO', 'BioRxiv'), ('CRO', 'Cross Ref')], max_length=3),
        ),
    ]
