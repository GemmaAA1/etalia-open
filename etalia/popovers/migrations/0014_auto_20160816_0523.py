# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('popovers', '0013_auto_20160816_0249'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='userperiodicemail',
            name='user',
        ),
        migrations.DeleteModel(
            name='UserPeriodicEmail',
        ),
    ]
