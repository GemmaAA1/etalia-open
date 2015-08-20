# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('feeds', '0019_userfeed_papers_seed'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfeed',
            name='message',
            field=models.CharField(blank=True, default='', max_length=127),
        ),
    ]
