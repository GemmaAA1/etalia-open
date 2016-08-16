# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0018_user_timezone'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userperiodicemail',
            name='last_sent_on',
            field=models.DateTimeField(null=True),
        ),
    ]
