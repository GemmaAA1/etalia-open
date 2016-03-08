# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0021_usersettings_stream_reactivity'),
    ]

    operations = [
        migrations.AddField(
            model_name='userlib',
            name='d_oldest',
            field=models.DateField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_author_weight',
            field=models.FloatField(default=0.1, verbose_name='Author weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_journal_weight',
            field=models.FloatField(default=0.1, verbose_name='Journal weight'),
        ),
        migrations.AlterField(
            model_name='usersettings',
            name='stream_reactivity',
            field=models.FloatField(default=0.5, verbose_name='Reactivity (higher reactivity focuses your stream on recent papers added to your library)'),
        ),
    ]
