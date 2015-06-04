# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='paper',
            name='type',
            field=models.CharField(default='', choices=[('JOU', 'Journal Article'), ('LET', 'Letter'), ('EDI', 'Editorial'), ('NEW', 'News'), ('PRO', 'Proceedings'), ('REV', 'Review'), ('PRE', 'e-Print'), ('DRA', 'Draft'), ('BOO', 'Book'), ('BOS', 'Book section'), ('PAT', 'Patent'), ('THE', 'Thesis'), ('', 'Unknown')], blank=True, max_length=3),
        ),
    ]
