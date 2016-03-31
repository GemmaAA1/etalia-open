# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from rest_framework.decorators import detail_route
from rest_framework.response import Response

from rest_condition import And, Or, Not
from etalia.core.api.permissions import IsReadOnlyRequest, IsPostRequest, \
    IsDeleteRequest, IsPutPatchRequest, IsThreadMember, IsOwner

from ..models import ThreadUser
from .serializers import ThreadUserSerializer

