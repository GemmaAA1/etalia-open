# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0024_auto_20150731_2313'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='FieldUseInModel',
            new_name='TextField',
        ),
        migrations.RenameField(
            model_name='textfield',
            old_name='field',
            new_name='text_field',
        ),
        migrations.RemoveField(
            model_name='model',
            name='fields',
        ),
        migrations.AddField(
            model_name='model',
            name='text_fields',
            field=models.ManyToManyField(related_name='text_fields', to='nlp.TextField'),
        ),
    ]
