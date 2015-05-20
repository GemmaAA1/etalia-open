# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import model_utils.fields
import core.models
import library.validators


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('last_name', models.CharField(default='', max_length=100, validators=[library.validators.validate_author_names])),
                ('first_name', models.CharField(default='', blank=True, max_length=100)),
                ('email', models.EmailField(default='', blank=True, max_length=254)),
            ],
        ),
        migrations.CreateModel(
            name='AuthorPaper',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='', blank=True, max_length=128, unique=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='CorpAuthorPaper',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('id_issn', core.models.NullableCharField(verbose_name='ISSN', max_length=9, unique=True, null=True, blank=True, default=None, validators=[library.validators.validate_issn])),
                ('id_eissn', core.models.NullableCharField(verbose_name='e-ISSN', max_length=9, unique=True, null=True, blank=True, default=None, validators=[library.validators.validate_issn])),
                ('id_arx', core.models.NullableCharField(verbose_name='Arxiv ID', max_length=32, unique=True, null=True, blank=True, default=None)),
                ('id_oth', core.models.NullableCharField(verbose_name='Other ID', max_length=32, unique=True, null=True, blank=True, default=None)),
                ('title', models.CharField(default='', max_length=200)),
                ('short_title', models.CharField(default='', blank=True, max_length=100)),
                ('url', models.URLField(default='', blank=True)),
                ('scope', models.TextField(default='', blank=True, max_length=1000)),
                ('language', models.CharField(default='ENG', blank=True, max_length=3, choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')])),
                ('period', models.CharField(default='', blank=True, max_length=200, choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular'), ('', 'Unknown')])),
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
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('publish_status', models.CharField(default='', blank=True, max_length=20, choices=[('ppublish', 'Paper Print'), ('epublish', 'e-Print only'), ('aheadofprint', 'e-Print ahead'), ('preprint', 'pre-Print'), ('', 'Unknown')])),
                ('publish_status_changed', model_utils.fields.MonitorField(monitor='publish_status', default=django.utils.timezone.now)),
                ('type', models.CharField(default='', blank=True, max_length=4, choices=[('JOUR', 'Journal Article'), ('LETT', 'Letter'), ('EDIT', 'Editorial'), ('NEWS', 'News'), ('PROC', 'Proceedings'), ('REVI', 'Review'), ('PREP', 'e-Print'), ('DRAF', 'Draft'), ('', 'Unknown')])),
                ('id_doi', core.models.NullableCharField(verbose_name='DOI', max_length=32, unique=True, null=True, blank=True, default='')),
                ('id_arx', core.models.NullableCharField(verbose_name='Arxiv', max_length=32, unique=True, null=True, blank=True, default='')),
                ('id_pmi', core.models.NullableCharField(verbose_name='PMID', max_length=32, unique=True, null=True, blank=True, default='')),
                ('id_pii', core.models.NullableCharField(verbose_name='PII', null=True, blank=True, max_length=32, default='')),
                ('id_oth', core.models.NullableCharField(verbose_name='Other ID', max_length=32, unique=True, null=True, blank=True, default='')),
                ('title', models.CharField(default='', max_length=500)),
                ('abstract', models.TextField(blank=True, default='')),
                ('volume', models.CharField(default='', blank=True, max_length=200)),
                ('issue', models.CharField(default='', blank=True, max_length=200)),
                ('page', models.CharField(default='', blank=True, max_length=200)),
                ('key_terms', models.TextField(blank=True, default='')),
                ('date_ep', models.DateField(null=True, blank=True, default=None)),
                ('date_pp', models.DateField(null=True, blank=True, default=None)),
                ('date_lr', models.DateField(null=True, blank=True, default=None)),
                ('url', models.URLField(default='', blank=True)),
                ('language', models.CharField(default='ENG', blank=True, max_length=3, choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')])),
                ('source', models.CharField(default='', blank=True, max_length=20)),
                ('source_changed', model_utils.fields.MonitorField(monitor='source', default=django.utils.timezone.now)),
                ('is_trusted', models.BooleanField(default=False)),
                ('authors', models.ManyToManyField(to='library.Author', blank=True, default=None, through='library.AuthorPaper')),
                ('corp_author', models.ManyToManyField(to='library.CorpAuthor', blank=True, default=None, through='library.CorpAuthorPaper')),
                ('journal', models.ForeignKey(to='library.Journal', null=True, blank=True, default=None)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(max_length=200, unique=True)),
                ('url', models.URLField(default=True, blank=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='journal',
            name='publisher',
            field=models.ForeignKey(to='library.Publisher', null=True, blank=True, default=None),
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
