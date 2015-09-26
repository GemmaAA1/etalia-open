# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0006_auto_20150926_0645'),
    ]

    operations = [
        migrations.AlterField(
            model_name='model',
            name='upload_state',
            field=models.CharField(max_length=3, default='IDL', choices=[('IDL', 'Idle'), ('ING', 'Uploading')]),
        ),
        migrations.AlterField(
            model_name='mostsimilar',
            name='upload_state',
            field=models.CharField(max_length=3, default='IDL ', choices=[('IDL', 'Idle'), ('ING', 'Uploading')]),
        ),
    ]
