# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0005_auto_20160625_0614'),
    ]

    operations = [
        migrations.AlterField(
            model_name='user',
            name='init_step',
            field=models.CharField(choices=[('NON', 'uninitialized'), ('LIB', 'library'), ('STR', 'paper feed'), ('TRE', 'trend feed'), ('THR', 'thread feed'), ('POP', 'popovers'), ('IDL', 'done')], help_text='Tag where init user stands', default='NON', max_length=3),
        ),
    ]
