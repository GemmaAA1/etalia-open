# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.contrib.postgres.fields


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0008_userfeedvector'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userfeedvector',
            name='feed',
            field=models.ForeignKey(to='feeds.UserFeed', related_name='vector'),
        ),
        migrations.AlterField(
            model_name='userfeedvector',
            name='vector',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.FloatField(null=True), size=300, null=True),
        ),
        migrations.AlterUniqueTogether(
            name='userfeedvector',
            unique_together=set([('feed', 'model')]),
        ),
    ]
