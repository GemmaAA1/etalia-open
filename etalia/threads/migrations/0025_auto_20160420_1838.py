# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('threads', '0024_auto_20160418_1841'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='thread',
            options={'ordering': ('-published_at',)},
        ),
        migrations.AddField(
            model_name='thread',
            name='published_at',
            field=models.DateTimeField(blank=True, null=True),
        ),
    ]
