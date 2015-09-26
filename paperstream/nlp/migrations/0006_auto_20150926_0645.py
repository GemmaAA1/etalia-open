# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0005_auto_20150925_0814'),
    ]

    operations = [
        migrations.AddField(
            model_name='model',
            name='upload_state',
            field=models.CharField(choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3, default='DON'),
        ),
        migrations.AddField(
            model_name='mostsimilar',
            name='upload_state',
            field=models.CharField(choices=[('IDL', 'Idle'), ('ING', 'Uploading')], max_length=3, default='DON'),
        ),
    ]
