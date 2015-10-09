# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0005_auto_20151009_0617'),
    ]

    operations = [
        migrations.RenameField(
            model_name='stats',
            old_name='nb_journal',
            new_name='nb_journals',
        ),
        migrations.RenameField(
            model_name='stats',
            old_name='nb_papers_last_two_week',
            new_name='nb_papers_last_two_weeks',
        ),
    ]
