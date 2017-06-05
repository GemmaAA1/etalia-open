# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0009_paper_date_co'),
    ]

    operations = [
        migrations.AddField(
            model_name='journal',
            name='is_in_fixture',
            field=models.BooleanField(default=True),
        ),
    ]
