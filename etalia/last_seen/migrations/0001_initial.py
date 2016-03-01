# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('sites', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='LastSeen',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', auto_created=True, serialize=False)),
                ('module', models.CharField(default='default', max_length=20)),
                ('last_seen', models.DateTimeField(default=django.utils.timezone.now)),
                ('site', models.ForeignKey(to='sites.Site')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ('-last_seen',),
            },
        ),
        migrations.AlterUniqueTogether(
            name='lastseen',
            unique_together=set([('user', 'site', 'module')]),
        ),
    ]
