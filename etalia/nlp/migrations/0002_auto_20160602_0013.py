# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('nlp', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('threads', '0001_initial'),
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='userfingerprint',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL, related_name='fingerprint'),
        ),
        migrations.AddField(
            model_name='threadvectors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='threadvectors',
            name='thread',
            field=models.ForeignKey(to='threads.Thread', related_name='vectors'),
        ),
        migrations.AddField(
            model_name='threadneighbors',
            name='te',
            field=models.ForeignKey(to='nlp.ThreadEngine'),
        ),
        migrations.AddField(
            model_name='threadneighbors',
            name='thread',
            field=models.ForeignKey(to='threads.Thread'),
        ),
        migrations.AddField(
            model_name='threadengine',
            name='model',
            field=models.ForeignKey(to='nlp.Model', related_name='te'),
        ),
        migrations.AddField(
            model_name='papervectors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='papervectors',
            name='paper',
            field=models.ForeignKey(to='library.Paper', related_name='vectors'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='paper',
            field=models.ForeignKey(to='library.Paper', related_name='neighbors'),
        ),
        migrations.AddField(
            model_name='paperneighbors',
            name='pe',
            field=models.ForeignKey(to='nlp.PaperEngine'),
        ),
        migrations.AddField(
            model_name='paperengine',
            name='model',
            field=models.ForeignKey(to='nlp.Model', related_name='pe'),
        ),
        migrations.AddField(
            model_name='model',
            name='text_fields',
            field=models.ManyToManyField(to='nlp.TextField', related_name='text_fields'),
        ),
        migrations.AddField(
            model_name='journalvectors',
            name='journal',
            field=models.ForeignKey(to='library.Journal', related_name='vectors'),
        ),
        migrations.AddField(
            model_name='journalvectors',
            name='model',
            field=models.ForeignKey(to='nlp.Model'),
        ),
        migrations.AddField(
            model_name='journalneighbors',
            name='journal',
            field=models.ForeignKey(to='library.Journal', related_name='neighbors'),
        ),
        migrations.AddField(
            model_name='journalneighbors',
            name='pe',
            field=models.ForeignKey(to='nlp.PaperEngine'),
        ),
        migrations.AlterUniqueTogether(
            name='userfingerprint',
            unique_together=set([('user', 'name')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadvectors',
            unique_together=set([('thread', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='threadneighbors',
            unique_together=set([('time_lapse', 'thread', 'te')]),
        ),
        migrations.AlterUniqueTogether(
            name='papervectors',
            unique_together=set([('paper', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='paperneighbors',
            unique_together=set([('time_lapse', 'paper', 'pe')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalvectors',
            unique_together=set([('journal', 'model')]),
        ),
        migrations.AlterUniqueTogether(
            name='journalneighbors',
            unique_together=set([('journal', 'pe')]),
        ),
    ]
