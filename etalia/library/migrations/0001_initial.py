# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import etalia.library.validators
import model_utils.fields
import etalia.core.models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('last_name', models.CharField(default='', validators=[etalia.library.validators.validate_author_names], max_length=100)),
                ('first_name', models.CharField(blank=True, default='', max_length=100)),
                ('email', models.EmailField(blank=True, default='', max_length=254)),
            ],
        ),
        migrations.CreateModel(
            name='AuthorPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('position', models.IntegerField(default=1)),
                ('author', models.ForeignKey(to='library.Author')),
            ],
            options={
                'ordering': ['position'],
            },
        ),
        migrations.CreateModel(
            name='CorpAuthor',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(blank=True, default='', unique=True, max_length=128)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='CorpAuthorPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('corp_author', models.ForeignKey(to='library.CorpAuthor')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Journal',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('id_issn', etalia.core.models.NullableCharField(blank=True, default='', unique=True, validators=[etalia.library.validators.validate_issn], max_length=9, db_index=True, verbose_name='ISSN', null=True)),
                ('id_eissn', etalia.core.models.NullableCharField(blank=True, default='', unique=True, validators=[etalia.library.validators.validate_issn], max_length=9, db_index=True, verbose_name='e-ISSN', null=True)),
                ('id_arx', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=32, null=True, verbose_name='Arxiv ID', db_index=True)),
                ('id_oth', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=32, null=True, verbose_name='Other ID', db_index=True)),
                ('title', models.CharField(default='', max_length=200)),
                ('short_title', models.CharField(blank=True, default='', max_length=100)),
                ('url', models.URLField(blank=True, default='')),
                ('scope', models.TextField(blank=True, default='', max_length=1000)),
                ('language', models.CharField(blank=True, default='ENG', choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')], max_length=3)),
                ('period', models.CharField(blank=True, default='', choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular'), ('', 'Unknown')], max_length=200)),
                ('is_trusted', models.BooleanField(default=False)),
                ('lib_size', models.IntegerField(default=0)),
            ],
            options={
                'ordering': ['title'],
            },
        ),
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('publish_status', models.CharField(blank=True, default='', choices=[('ppublish', 'Paper Print'), ('epublish', 'e-Print only'), ('aheadofprint', 'e-Print ahead'), ('preprint', 'pre-Print'), ('', 'Unknown')], max_length=20)),
                ('publish_status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='publish_status')),
                ('type', models.CharField(blank=True, default='', choices=[('JOU', 'Journal Article'), ('LET', 'Letter'), ('EDI', 'Editorial'), ('NEW', 'News'), ('PRO', 'Proceedings'), ('REV', 'Review'), ('PRE', 'e-Print'), ('DRA', 'Draft'), ('BOO', 'Book'), ('BOS', 'Book section'), ('PAT', 'Patent'), ('THE', 'Thesis'), ('', 'Unknown')], max_length=3)),
                ('id_doi', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=64, null=True, verbose_name='DOI', db_index=True)),
                ('id_arx', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=64, null=True, verbose_name='Arxiv', db_index=True)),
                ('id_pmi', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=64, null=True, verbose_name='PMID', db_index=True)),
                ('id_pii', etalia.core.models.NullableCharField(blank=True, default='', max_length=64, null=True, verbose_name='PII', db_index=True)),
                ('id_isbn', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=64, null=True, verbose_name='ISBN', db_index=True)),
                ('id_oth', etalia.core.models.NullableCharField(blank=True, default='', unique=True, max_length=64, null=True, verbose_name='Other ID', db_index=True)),
                ('title', models.CharField(default='', max_length=500, db_index=True)),
                ('abstract', models.TextField(blank=True, default='')),
                ('volume', models.CharField(blank=True, default='', max_length=200)),
                ('issue', models.CharField(blank=True, default='', max_length=200)),
                ('page', models.CharField(blank=True, default='', max_length=200)),
                ('key_terms', models.TextField(blank=True, default='')),
                ('date_ep', models.DateField(blank=True, default=None, db_index=True, null=True)),
                ('date_pp', models.DateField(blank=True, default=None, db_index=True, null=True)),
                ('date_lr', models.DateField(blank=True, default=None, null=True)),
                ('date_fs', models.DateField(auto_now_add=True, db_index=True)),
                ('url', models.URLField(blank=True, default='')),
                ('language', models.CharField(blank=True, default='ENG', choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')], max_length=3)),
                ('source', models.CharField(blank=True, default='', max_length=20)),
                ('source_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='source')),
                ('is_trusted', models.BooleanField(default=False, db_index=True)),
                ('authors', models.ManyToManyField(blank=True, to='library.Author', default=None, through='library.AuthorPaper')),
                ('corp_author', models.ManyToManyField(blank=True, to='library.CorpAuthor', default=None, through='library.CorpAuthorPaper')),
                ('journal', models.ForeignKey(blank=True, default=None, null=True, to='library.Journal')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(unique=True, max_length=200)),
                ('url', models.URLField(blank=True, default=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Stats',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('nb_papers', models.IntegerField(default=0)),
                ('nb_journals', models.IntegerField(default=0)),
                ('nb_authors', models.IntegerField(default=0)),
                ('nb_papers_last_week', models.IntegerField(default=0)),
                ('nb_papers_last_two_weeks', models.IntegerField(default=0)),
                ('nb_papers_last_month', models.IntegerField(default=0)),
                ('nb_papers_last_year', models.IntegerField(default=0)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(blank=True, default=None, null=True, to='library.Publisher'),
        ),
        migrations.AddField(
            model_name='corpauthorpaper',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AddField(
            model_name='authorpaper',
            name='paper',
            field=models.ForeignKey(to='library.Paper'),
        ),
        migrations.AlterUniqueTogether(
            name='author',
            unique_together=set([('first_name', 'last_name')]),
        ),
    ]
