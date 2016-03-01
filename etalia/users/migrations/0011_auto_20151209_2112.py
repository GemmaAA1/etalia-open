# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0003_journal_lib_size'),
        ('users', '0010_auto_20151209_0127'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserLibAuthor',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('occurrence', models.IntegerField(default=0)),
                ('author', models.ForeignKey(to='library.Author')),
                ('userlib', models.ForeignKey(to='users.UserLib')),
            ],
            options={
                'ordering': ('-occurrence',),
            },
        ),
        migrations.AddField(
            model_name='userlib',
            name='authors',
            field=models.ManyToManyField(to='library.Author', through='users.UserLibAuthor'),
        ),
        migrations.AlterUniqueTogether(
            name='userlibauthor',
            unique_together=set([('userlib', 'author')]),
        ),
    ]
