# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0005_auto_20161007_1913'),
    ]

    operations = [
        migrations.CreateModel(
            name='ConsumerSpringer',
            fields=[
                ('consumer_ptr', models.OneToOneField(serialize=False, parent_link=True, to='consumers.Consumer', auto_created=True, primary_key=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('consumers.consumer',),
        ),
        migrations.AlterField(
            model_name='consumer',
            name='type',
            field=models.CharField(max_length=3, choices=[('PUB', 'PubMed'), ('ELS', 'Elsevier'), ('ARX', 'Arxiv'), ('IEE', 'IEEE'), ('SPR', 'Springer'), ('BIO', 'BioRxiv'), ('CRO', 'Cross Ref')]),
        ),
    ]
