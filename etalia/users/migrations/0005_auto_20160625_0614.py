# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0004_auto_20160605_2039'),
    ]

    operations = [
        migrations.RenameField(
            model_name='usersettings',
            old_name='stream_roll_back_deltatime',
            new_name='fingerprint_roll_back_deltatime',
        ),
        migrations.AlterField(
            model_name='user',
            name='init_step',
            field=models.CharField(default='NON', max_length=3, choices=[('NON', 'uninitialized'), ('LIB', 'library'), ('STR', 'paper feed'), ('TRE', 'trend feed'), ('THR', 'thread feed'), ('IDL', 'done')], help_text='Tag where init user stands'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_vector_weight',
            field=models.FloatField(default=1.0, verbose_name='Title/Abstract weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_altmetric_weight',
            field=models.FloatField(default=0.5, verbose_name='Altmetric weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='trend_doc_weight',
            field=models.FloatField(default=1.0, verbose_name='Title/Abstract weight'),
        ),
    ]
