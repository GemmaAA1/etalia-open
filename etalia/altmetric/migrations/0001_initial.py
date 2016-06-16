# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('library', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='AltmetricModel',
            fields=[
                ('id', models.AutoField(primary_key=True, verbose_name='ID', serialize=False, auto_created=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('score', models.FloatField(default=0.0, db_index=True)),
                ('altmetric_id', models.IntegerField(default=0)),
                ('altmetric_jid', models.CharField(blank=True, max_length=32)),
                ('score_1d', models.FloatField(default=0.0)),
                ('score_2d', models.FloatField(default=0.0)),
                ('score_3d', models.FloatField(default=0.0)),
                ('score_4d', models.FloatField(default=0.0)),
                ('score_5d', models.FloatField(default=0.0)),
                ('score_6d', models.FloatField(default=0.0)),
                ('score_1w', models.FloatField(default=0.0)),
                ('score_1m', models.FloatField(default=0.0)),
                ('score_3m', models.FloatField(default=0.0)),
                ('score_6m', models.FloatField(default=0.0)),
                ('score_1y', models.FloatField(default=0.0)),
                ('cited_by_posts_count', models.IntegerField(default=0)),
                ('cited_by_delicious_count', models.IntegerField(default=0)),
                ('cited_by_fbwalls_count', models.IntegerField(default=0)),
                ('cited_by_feeds_count', models.IntegerField(default=0)),
                ('cited_by_forum_count', models.IntegerField(default=0)),
                ('cited_by_gplus_count', models.IntegerField(default=0)),
                ('cited_by_linkedin_count', models.IntegerField(default=0)),
                ('cited_by_msm_count', models.IntegerField(default=0)),
                ('cited_by_peer_review_sites_count', models.IntegerField(default=0)),
                ('cited_by_pinners_count', models.IntegerField(default=0)),
                ('cited_by_policies_count', models.IntegerField(default=0)),
                ('cited_by_qs_count', models.IntegerField(default=0)),
                ('cited_by_rdts_count', models.IntegerField(default=0)),
                ('cited_by_rh_count', models.IntegerField(default=0)),
                ('cited_by_tweeters_count', models.IntegerField(default=0)),
                ('cited_by_videos_count', models.IntegerField(default=0)),
                ('cited_by_weibo_count', models.IntegerField(default=0)),
                ('cited_by_wikipedia_count', models.IntegerField(default=0)),
                ('readers_citeulike', models.IntegerField(default=0)),
                ('readers_mendeley', models.IntegerField(default=0)),
                ('image', models.URLField(default='')),
                ('type', models.CharField(default='zzzzzzzz', max_length=8)),
                ('paper', models.OneToOneField(to='library.Paper', related_name='altmetric')),
            ],
            options={
                'ordering': ['-score'],
            },
        ),
    ]
