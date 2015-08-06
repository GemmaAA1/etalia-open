# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0004_paper_id_isbn'),
        ('nlp', '0005_auto_20150708_1632'),
    ]

    operations = [
        migrations.CreateModel(
            name='PaperVectors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('vector', django.contrib.postgres.fields.ArrayField(size=None, base_field=models.FloatField())),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='fingerprint',
            name='mod',
        ),
        migrations.RemoveField(
            model_name='fingerprint',
            name='paper',
        ),
        migrations.AlterField(
            model_name='model',
            name='status',
            field=models.CharField(default='UNT', choices=[('UNT', 'Untrained'), ('VOC', 'Building Vocabulary'), ('TRA', 'Training'), ('SAV', 'Saving'), ('LOA', 'Loading'), ('IDL', 'Idle'), ('USE', 'Usable')], max_length=3),
        ),
        migrations.DeleteModel(
            name='Fingerprint',
        ),
        migrations.AddField(
            model_name='papervectors',
            name='mod',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='papervectors',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
    ]
