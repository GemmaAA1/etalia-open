from django.db import models
from django.conf import settings

from etalia.core.models import TimeStampedModel
from .constants import POPOVER_TYPES, POPOVER_STATUSES


class PopOver(TimeStampedModel):

    title = models.CharField(max_length=256)

    body = models.TextField()

    anchor = models.CharField(max_length=128)

    type = models.PositiveIntegerField(choices=POPOVER_TYPES, default=1)


class UserPopOver(TimeStampedModel):

    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    popover = models.ForeignKey(PopOver)

    status = models.PositiveIntegerField(choices=POPOVER_STATUSES,
                                         default=1)

