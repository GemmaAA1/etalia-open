# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0008_auto_20150506_0204'),
    ]

    operations = [
        migrations.CreateModel(
            name='AuthorPosition',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('position', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'ordering': ['position'],
            },
        ),
        migrations.RenameModel(
            old_name='Creator',
            new_name='Author',
        ),
        migrations.RemoveField(
            model_name='creatorposition',
            name='creator',
        ),
        migrations.RemoveField(
            model_name='creatorposition',
            name='paper',
        ),
        migrations.RemoveField(
            model_name='paper',
            name='creators',
        ),
        migrations.DeleteModel(
            name='CreatorPosition',
        ),
        migrations.AddField(
            model_name='authorposition',
            name='author',
            field=models.ForeignKey(to='library.Author'),
        ),
        migrations.AddField(
            model_name='authorposition',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AddField(
            model_name='paper',
            name='authors',
            field=models.ManyToManyField(through='library.AuthorPosition', to='library.Author'),
        ),
    ]
