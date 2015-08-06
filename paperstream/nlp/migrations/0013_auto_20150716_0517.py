# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0012_auto_20150711_0157'),
    ]

    operations = [
        migrations.CreateModel(
            name='FieldUseInModel',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('field', models.CharField(choices=[('title', 'Title'), ('abstract', 'Abstract')], max_length=100)),
            ],
        ),
        migrations.AddField(
            model_name='model',
            name='data_path',
            field=models.CharField(default='/Users/nicolaspannetier/Projects/paperstream/paperstream_project/paperstream/nlp/data', max_length=256),
        ),
        migrations.AddField(
            model_name='fielduseinmodel',
            name='model',
            field=models.ForeignKey(related_name='paperfields', to='nlp.Model'),
        ),
        migrations.AlterUniqueTogether(
            name='fielduseinmodel',
            unique_together=set([('model', 'field')]),
        ),
    ]
