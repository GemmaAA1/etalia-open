# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0020_auto_20161104_1812'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='init_step',
            field=models.CharField(default='NON', help_text='Tag where init user stands', max_length=3, choices=[('NON', 'uninitialized'), ('LIB', 'Syncing Papers'), ('STR', 'Syncing Feed'), ('TRE', 'Syncing Trend'), ('THR', 'Syncing Thread'), ('POP', 'Finalizing'), ('IDL', 'done')]),
        ),
    ]
