# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0032_auto_20150515_0022'),
    ]

    operations = [
        migrations.AlterField(
            model_name='journal',
            name='short_title',
            field=models.CharField(max_length=100, blank=True, default=''),
        ),
        migrations.AlterField(
            model_name='journal',
            name='title',
            field=models.CharField(max_length=200, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='authors',
            field=models.ManyToManyField(through='library.AuthorPaper', to='library.Author', blank=True, default=None),
        ),
        migrations.AlterField(
            model_name='paper',
            name='corp_author',
            field=models.ManyToManyField(through='library.CorpAuthorPaper', to='library.CorpAuthor', blank=True, default=None),
        ),
        migrations.AlterField(
            model_name='paper',
            name='date_ep',
            field=models.DateField(null=True, blank=True, default=None),
        ),
        migrations.AlterField(
            model_name='paper',
            name='date_lr',
            field=models.DateField(null=True, blank=True, default=None),
        ),
        migrations.AlterField(
            model_name='paper',
            name='date_pp',
            field=models.DateField(null=True, blank=True, default=None),
        ),
        migrations.AlterField(
            model_name='paper',
            name='journal',
            field=models.ForeignKey(to='library.Journal', null=True, blank=True, default=None),
        ),
        migrations.AlterField(
            model_name='paper',
            name='source',
            field=models.CharField(max_length=20, blank=True, choices=[('MEND', 'Mendeley'), ('ARXI', 'Arxiv'), ('ZOTE', 'Zotero'), ('ELSV', 'Elsevier'), ('PUMD', 'PubMed'), ('IEEE', 'IEEE'), ('', 'Unknown')], default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='title',
            field=models.CharField(max_length=500, default=''),
        ),
    ]
