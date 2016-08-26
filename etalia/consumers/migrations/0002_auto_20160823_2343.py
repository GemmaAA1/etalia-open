# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PubPeerConsumer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AlterField(
            model_name='consumer',
            name='type',
            field=models.CharField(max_length=3, choices=[('PUB', 'PubMed'), ('ELS', 'Elsevier'), ('ARX', 'Arxiv'), ('IEE', 'IEEE'), ('CRO', 'Cross Ref')]),
        ),
    ]
