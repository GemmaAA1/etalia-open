# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models
import model_utils.fields
import django.utils.timezone
import library.validators


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Author',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('last_name', models.CharField(validators=[library.validators.validate_author_names], default='', max_length=100)),
                ('first_name', models.CharField(default='', max_length=100, blank=True)),
                ('email', models.EmailField(default='', max_length=254, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='AuthorPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
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
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(default='', max_length=128, unique=True, blank=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='CorpAuthorPaper',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
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
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('id_issn', core.models.NullableCharField(validators=[library.validators.validate_issn], null=True, blank=True, verbose_name='ISSN', default=None, max_length=9, unique=True)),
                ('id_eissn', core.models.NullableCharField(validators=[library.validators.validate_issn], null=True, blank=True, verbose_name='e-ISSN', default=None, max_length=9, unique=True)),
                ('id_arx', core.models.NullableCharField(null=True, blank=True, verbose_name='Arxiv ID', default=None, max_length=32, unique=True)),
                ('id_oth', core.models.NullableCharField(null=True, blank=True, verbose_name='Other ID', default=None, max_length=32, unique=True)),
                ('title', models.CharField(default='', max_length=200)),
                ('short_title', models.CharField(default='', max_length=100, blank=True)),
                ('url', models.URLField(default='', blank=True)),
                ('scope', models.TextField(default='', max_length=1000, blank=True)),
                ('language', models.CharField(choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')], default='ENG', max_length=3, blank=True)),
                ('period', models.CharField(choices=[('ANN', 'Annual'), ('SEM', 'Semi-annual'), ('TRI', 'Tri-annual'), ('QUA', 'Quarterly'), ('MON', 'Monthly'), ('BIM', 'Bi-monthly'), ('IRR', 'Irregular'), ('', 'Unknown')], default='', max_length=200, blank=True)),
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
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('publish_status', models.CharField(choices=[('ppublish', 'Paper Print'), ('epublish', 'e-Print only'), ('aheadofprint', 'e-Print ahead'), ('preprint', 'pre-Print'), ('', 'Unknown')], default='', max_length=20, blank=True)),
                ('publish_status_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='publish_status')),
                ('type', models.CharField(choices=[('JOU', 'Journal Article'), ('LET', 'Letter'), ('EDI', 'Editorial'), ('NEW', 'News'), ('PRO', 'Proceedings'), ('REV', 'Review'), ('PRE', 'e-Print'), ('DRA', 'Draft'), ('', 'Unknown')], default='', max_length=3, blank=True)),
                ('id_doi', core.models.NullableCharField(null=True, blank=True, verbose_name='DOI', default='', max_length=64, unique=True)),
                ('id_arx', core.models.NullableCharField(null=True, blank=True, verbose_name='Arxiv', default='', max_length=64, unique=True)),
                ('id_pmi', core.models.NullableCharField(null=True, blank=True, verbose_name='PMID', default='', max_length=64, unique=True)),
                ('id_pii', core.models.NullableCharField(verbose_name='PII', default='', max_length=64, null=True, blank=True)),
                ('id_oth', core.models.NullableCharField(null=True, blank=True, verbose_name='Other ID', default='', max_length=64, unique=True)),
                ('title', models.CharField(default='', max_length=500)),
                ('abstract', models.TextField(default='', blank=True)),
                ('volume', models.CharField(default='', max_length=200, blank=True)),
                ('issue', models.CharField(default='', max_length=200, blank=True)),
                ('page', models.CharField(default='', max_length=200, blank=True)),
                ('key_terms', models.TextField(default='', blank=True)),
                ('date_ep', models.DateField(default=None, null=True, blank=True)),
                ('date_pp', models.DateField(default=None, null=True, blank=True)),
                ('date_lr', models.DateField(default=None, null=True, blank=True)),
                ('url', models.URLField(default='', blank=True)),
                ('language', models.CharField(choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh')], default='ENG', max_length=3, blank=True)),
                ('source', models.CharField(default='', max_length=20, blank=True)),
                ('source_changed', model_utils.fields.MonitorField(default=django.utils.timezone.now, monitor='source')),
                ('is_trusted', models.BooleanField(default=False)),
                ('authors', models.ManyToManyField(through='library.AuthorPaper', default=None, to='library.Author', blank=True)),
                ('corp_author', models.ManyToManyField(through='library.CorpAuthorPaper', default=None, to='library.CorpAuthor', blank=True)),
                ('journal', models.ForeignKey(null=True, blank=True, default=None, to='library.Journal')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Publisher',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, primary_key=True, auto_created=True)),
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
            field=models.ForeignKey(null=True, blank=True, default=None, to='library.Publisher'),
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
