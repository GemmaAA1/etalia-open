# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.db.models import Prefetch

from etalia.users.constants import USER_INDIVIDUAL

from .models import ThreadUser


class ThreadMixin(object):

    @staticmethod
    def setup_eager_loading(queryset, **kwargs):
        """ Perform necessary eager loading of data. """
        queryset = queryset.select_related('user')
        if 'user'in kwargs and \
            kwargs['user'].is_authenticated() and \
            kwargs['user'].type == USER_INDIVIDUAL:
            queryset = queryset.prefetch_related(
                Prefetch('threaduser_set',
                         to_attr='tu',
                         queryset=ThreadUser.objects.filter(user=kwargs['user'])))
        return queryset