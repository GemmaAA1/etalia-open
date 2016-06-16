# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import etalia.core.mixins
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PaperUser',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, auto_created=True, verbose_name='ID')),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('watch', models.PositiveIntegerField(default=None, null=True, choices=[(1, 'Pinned'), (2, 'Banned')])),
            ],
            bases=(etalia.core.mixins.ModelDiffMixin, models.Model),
        ),
        migrations.CreateModel(
            name='PaperUserHistory',
            fields=[
                ('paperuser_ptr', models.OneToOneField(serialize=False, to='library.PaperUser', primary_key=True, auto_created=True, parent_link=True)),
                ('difference', models.CharField(default='', max_length=256)),
                ('date', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'abstract': False,
            },
            bases=('library.paperuser',),
        ),
        migrations.AddField(
            model_name='paperuser',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AddField(
            model_name='paperuser',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='paperuserhistory',
            name='paperuser',
            field=models.ForeignKey(to='library.PaperUser', related_name='history'),
        ),
        migrations.AlterUniqueTogether(
            name='paperuser',
            unique_together=set([('paper', 'user')]),
        ),
    ]
