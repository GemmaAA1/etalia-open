# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework import serializers
from django.contrib.auth import get_user_model

User = get_user_model()


class UserSerializer(serializers.ModelSerializer):
    url = serializers.URLField(source='get_absolute_url')

    class Meta:
        model = User
        fields = ('pk',
                  'email',
                  'first_name',
                  'last_name',
                  'url')
