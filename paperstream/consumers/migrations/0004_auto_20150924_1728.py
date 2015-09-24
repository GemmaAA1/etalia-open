# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('consumers', '0003_auto_20150914_0029'),
    ]

    operations = [
        migrations.AddField(
            model_name='consumerjournalstat',
            name='message',
            field=models.CharField(default='', max_length=512),
        ),
        migrations.AlterField(
            model_name='consumerjournalstat',
            name='number_papers_fetched',
            field=models.IntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='consumerjournalstat',
            name='number_papers_recorded',
            field=models.IntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='consumerjournalstat',
            name='status',
            field=models.CharField(max_length=3, choices=[('SUC', 'Success'), ('FAI', 'Failed'), ('RES', 'Reset'), ('ACT', 'Activate'), ('DEA', 'Deactivate')]),
        ),
    ]
