# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.shortcuts import resolve_url
from rest_framework import exceptions, status


class MendeleyRedirectLoginErrorSerializer(exceptions.APIException):
    status_code = status.HTTP_302_FOUND

    def __init__(self):
        self.detail = resolve_url('social:begin', backend='mendeley')

    def __str__(self):
        return self.detail
