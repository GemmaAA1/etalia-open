# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import model_utils.fields
import core.models
import django.utils.timezone
import library.validators


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('last_name', models.CharField(max_length=100, default='', validators=[library.validators.validate_author_names])),
                ('first_name', models.CharField(blank=True, max_length=100, default='')),
                ('email', models.EmailField(blank=True, max_length=254, default='')),
            ],
        ),
        migrations.CreateModel(
            name='AuthorPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(blank=True, max_length=128, default='', unique=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='CorpAuthorPaper',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
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
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('id_issn', core.models.NullableCharField(blank=True, max_length=9, verbose_name='ISSN', db_index=True, unique=True, null=True, default='', validators=[library.validators.validate_issn])),
                ('id_eissn', core.models.NullableCharField(blank=True, max_length=9, verbose_name='e-ISSN', db_index=True, unique=True, null=True, default='', validators=[library.validators.validate_issn])),
                ('id_arx', core.models.NullableCharField(blank=True, max_length=32, verbose_name='Arxiv ID', db_index=True, unique=True, null=True, default='')),
                ('id_oth', core.models.NullableCharField(blank=True, max_length=32, verbose_name='Other ID', db_index=True, unique=True, null=True, default='')),
                ('title', models.CharField(max_length=200, default='')),
                ('short_title', models.CharField(blank=True, max_length=100, default='')),
                ('url', models.URLField(blank=True, default='')),
                ('scope', models.TextField(blank=True, max_length=1000, default='')),
                ('language', models.CharField(choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')], max_length=3, default='ENG', blank=True)),
                ('period', models.CharField(choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular'), ('', 'Unknown')], max_length=200, default='', blank=True)),
                ('lib_size', models.IntegerField(default=0)),
                ('is_trusted', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['title'],
            },
        ),
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('publish_status', models.CharField(choices=[('ppublish', 'Paper Print'), ('epublish', 'e-Print only'), ('aheadofprint', 'e-Print ahead'), ('preprint', 'pre-Print'), ('', 'Unknown')], max_length=20, default='', blank=True)),
                ('publish_status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='publish_status')),
                ('type', models.CharField(choices=[('JOU', 'Journal Article'), ('LET', 'Letter'), ('EDI', 'Editorial'), ('NEW', 'News'), ('PRO', 'Proceedings'), ('REV', 'Review'), ('PRE', 'e-Print'), ('DRA', 'Draft'), ('BOO', 'Book'), ('BOS', 'Book section'), ('PAT', 'Patent'), ('THE', 'Thesis'), ('', 'Unknown')], max_length=3, default='', blank=True)),
                ('id_doi', core.models.NullableCharField(blank=True, max_length=64, verbose_name='DOI', db_index=True, unique=True, null=True, default='')),
                ('id_arx', core.models.NullableCharField(blank=True, max_length=64, verbose_name='Arxiv', db_index=True, unique=True, null=True, default='')),
                ('id_pmi', core.models.NullableCharField(blank=True, max_length=64, verbose_name='PMID', db_index=True, unique=True, null=True, default='')),
                ('id_pii', core.models.NullableCharField(blank=True, max_length=64, verbose_name='PII', db_index=True, null=True, default='')),
                ('id_isbn', core.models.NullableCharField(blank=True, max_length=64, verbose_name='ISBN', db_index=True, unique=True, null=True, default='')),
                ('id_oth', core.models.NullableCharField(blank=True, max_length=64, verbose_name='Other ID', db_index=True, unique=True, null=True, default='')),
                ('title', models.CharField(max_length=500, default='', db_index=True)),
                ('abstract', models.TextField(blank=True, default='')),
                ('volume', models.CharField(blank=True, max_length=200, default='')),
                ('issue', models.CharField(blank=True, max_length=200, default='')),
                ('page', models.CharField(blank=True, max_length=200, default='')),
                ('key_terms', models.TextField(blank=True, default='')),
                ('date_ep', models.DateField(null=True, blank=True, default=None, db_index=True)),
                ('date_pp', models.DateField(null=True, blank=True, default=None, db_index=True)),
                ('date_lr', models.DateField(null=True, blank=True, default=None)),
                ('url', models.URLField(blank=True, default='')),
                ('language', models.CharField(choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')], max_length=3, default='ENG', blank=True)),
                ('source', models.CharField(blank=True, max_length=20, default='')),
                ('source_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='source')),
                ('is_trusted', models.BooleanField(default=False)),
                ('authors', models.ManyToManyField(to='library.Author', blank=True, through='library.AuthorPaper', default=None)),
                ('corp_author', models.ManyToManyField(to='library.CorpAuthor', blank=True, through='library.CorpAuthorPaper', default=None)),
                ('journal', models.ForeignKey(to='library.Journal', blank=True, null=True, default=None)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(primary_key=True, auto_created=True, verbose_name='ID', serialize=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(unique=True, max_length=200)),
                ('url', models.URLField(blank=True, default=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(to='library.Publisher', blank=True, null=True, default=None),
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
