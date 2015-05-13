# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0004_auto_20150513_0605'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='arxivconsumer',
            name='url',
        ),
        migrations.AlterUniqueTogether(
            name='consumerjournal',
            unique_together=set([('journal', 'consumer')]),
        ),
    ]
