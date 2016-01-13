# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0015_auto_20160112_0736'),
    ]

    operations = [
        migrations.AddField(
            model_name='usersettings',
            name='stream_author_weight',
            field=models.FloatField(default=1.0),
        ),
        migrations.AddField(
            model_name='usersettings',
            name='stream_journal_weight',
            field=models.FloatField(default=1.0),
        ),
        migrations.AddField(
            model_name='usersettings',
            name='stream_vector_weight',
            field=models.FloatField(default=1.0),
        ),
        migrations.AddField(
            model_name='usersettings',
            name='trend_altmetric_weight',
            field=models.FloatField(default=1.0),
        ),
        migrations.AddField(
            model_name='usersettings',
            name='trend_doc_weight',
            field=models.FloatField(default=1.0),
        ),
    ]
