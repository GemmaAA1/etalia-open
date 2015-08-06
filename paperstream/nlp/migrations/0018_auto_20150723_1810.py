# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0017_auto_20150723_1751'),
    ]

    operations = [
        migrations.CreateModel(
            name='LSH',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('state', models.CharField(max_length=3, default='', choices=[('NON', 'None'), ('BUS', 'Busy'), ('IDL', 'Idle')])),
                ('model', models.OneToOneField(to='nlp.Model', related_name='lsh')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='lshfmodel',
            name='model',
        ),
        migrations.RenameField(
            model_name='papervectors',
            old_name='is_in_lshf',
            new_name='is_in_full_lsh',
        ),
        migrations.DeleteModel(
            name='LSHFModel',
        ),
    ]
