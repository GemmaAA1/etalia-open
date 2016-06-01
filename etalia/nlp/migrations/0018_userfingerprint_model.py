# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0017_userfingerprint'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfingerprint',
            name='model',
            field=models.ForeignKey(default=1, to='nlp.Model'),
            preserve_default=False,
        ),
    ]
