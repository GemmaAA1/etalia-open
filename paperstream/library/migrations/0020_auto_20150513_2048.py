# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import core.models


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0019_auto_20150513_2013'),
    ]

    operations = [
        migrations.AddField(
            model_name='paper',
            name='source',
            field=models.CharField(blank=True, choices=[('PBMD', 'PubMed'), ('ELSV', 'Elsevier'), ('ARXI', 'Arxiv'), ('IEEE', 'IEEE'), ('ZOTE', 'Zotero'), ('MEND', 'Mendeley'), ('UNKN', 'Unknown')], max_length=4),
        ),
        migrations.AlterField(
            model_name='journal',
            name='language',
            field=models.CharField(blank=True, choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('PUS', 'Pushto'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh'), ('', 'Unknown')], default='ENG', max_length=3),
        ),
        migrations.AlterField(
            model_name='paper',
            name='id_pii',
            field=core.models.NullableCharField(blank=True, unique=True, verbose_name='PII', max_length=32, null=True, default=''),
        ),
        migrations.AlterField(
            model_name='paper',
            name='language',
            field=models.CharField(blank=True, choices=[('ENG', 'English'), ('AFR', 'Afrikaans'), ('ALB', 'Albanian'), ('AMH', 'Amharic'), ('ARA', 'Arabic'), ('ARM', 'Armenian'), ('AZE', 'Azerbaijani'), ('BEN', 'Bengali'), ('BOS', 'Bosnian'), ('BUL', 'Bulgarian'), ('CAT', 'Catalan'), ('CHI', 'Chinese'), ('CZE', 'Czech'), ('DAN', 'Danish'), ('DUT', 'Dutch'), ('ENG', 'English'), ('EPO', 'Esperanto'), ('EST', 'Estonian'), ('FIN', 'Finnish'), ('FRE', 'French'), ('GEO', 'Georgian'), ('GER', 'German'), ('GLA', 'Scottish Gaelic'), ('GRE', 'Greek, Modern'), ('HEB', 'Hebrew'), ('HIN', 'Hindi'), ('HRV', 'Croatian'), ('HUN', 'Hungarian'), ('ICE', 'Icelandic'), ('IND', 'Indonesian'), ('ITA', 'Italian'), ('JPN', 'Japanese'), ('KIN', 'Kinyarwanda'), ('KOR', 'Korean'), ('LAT', 'Latin'), ('LAV', 'Latvian'), ('LIT', 'Lithuanian'), ('MAC', 'Macedonian'), ('MAL', 'Malayalam'), ('MAO', 'Maori'), ('MAY', 'Malay'), ('MUL', 'Multiple languages'), ('NOR', 'Norwegian'), ('PER', 'Persian, Iranian'), ('POL', 'Polish'), ('POR', 'Portuguese'), ('PUS', 'Pushto'), ('RUM', 'Romanian, Rumanian, Moldovan'), ('RUS', 'Russian'), ('SAN', 'Sanskrit'), ('SLO', 'Slovak'), ('SLV', 'Slovenian'), ('SPA', 'Spanish'), ('SRP', 'Serbian'), ('SWE', 'Swedish'), ('THA', 'Thai'), ('TUR', 'Turkish'), ('UKR', 'Ukrainian'), ('UND', 'Undetermined'), ('URD', 'Urdu'), ('VIE', 'Vietnamese'), ('WEL', 'Welsh'), ('', 'Unknown')], default='ENG', max_length=3),
        ),
    ]
