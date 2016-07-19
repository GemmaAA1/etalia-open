
import json
import os
from django.db import models, transaction
from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.validators import MaxValueValidator, MinValueValidator

from etalia.core.models import TimeStampedModel
from .constants import POPOVER_TYPES, POPOVER_STATUSES, NEW, \
    POPOVER_DISPLAY_CHOICES, DISPLAY, MODAL, ANCHORED, HIDE, GOT_IT


class PopOver(TimeStampedModel):

    title = models.CharField(max_length=256, null=True, blank=True)

    template_path = models.CharField(max_length=128, null=True, blank=True)

    anchor = models.CharField(max_length=128, null=True, blank=True)

    type = models.PositiveIntegerField(choices=POPOVER_TYPES, default=1)

    # Priority: Priority scale from 1 (highest) to 9 (lowest)
    priority = models.PositiveIntegerField(default=1,
                                           validators=[MinValueValidator(1),
                                                       MaxValueValidator(9), ])

    def __str__(self):
        return self.title

    def reset(self):
        UserPopOver.object.filter(popover=self).update(status=NEW)

    def init(self):
        User = get_user_model()
        us = User.objects.all()
        objs = []
        UserPopOver.objects.filter(popover=self).delete()
        for user in us:
            objs.append(UserPopOver(user=user, popover=self, status=NEW))
        # bulk create UserPopOver
        UserPopOver.objects.bulk_create(objs)

    @classmethod
    def load_from_file(cls,
                       file=os.path.join(os.path.dirname(__file__), 'popovers.json')):

        with open(file) as data_file:
            data = json.load(data_file)['popovers']

        # Check that id field is unique
        ids = []
        for po in data:
            ids.append(po.get('id'))
        assert len(ids) == len(set(ids))

        # Reorder popovers
        new_pos = {po.pop('id'): po for po in data}

        # Update popovers
        for id_, new_po in new_pos.items():
            po, new = cls.objects.get_or_create(id=id_)
            for attr, val in new_po.items():
                setattr(po, attr, val)
            po.save()
            if new:   # init UserPopOver
                po.init()

        # Update display
        User = get_user_model()
        us = User.objects.all()
        for user in us:
            upoud, _ = UserPopOverUpdateDisplay.objects.get_or_create(user=user)
            upoud.update_display()


class UserPopOver(TimeStampedModel):

    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    popover = models.ForeignKey(PopOver)

    status = models.PositiveIntegerField(choices=POPOVER_STATUSES,
                                         default=NEW)

    display = models.PositiveIntegerField(choices=POPOVER_DISPLAY_CHOICES,
                                          default=HIDE)

    class Meta:
        unique_together = ('popover', 'user')
        ordering = ('-popover__type', 'popover__priority')

    def __str__(self):
        return self.popover.title


class UserPopOverUpdateDisplay(TimeStampedModel):

    user = models.ForeignKey(settings.AUTH_USER_MODEL)

    task_id = models.CharField(max_length=128, null=True, blank=True,
                               default='')

    def deferred_display_update(self):
        from .tasks import update_popovers_display

        # trigger update_display deferred task if not defined
        if not self.task_id:
            # plan display_update
            res = update_popovers_display.apply_async(
                args=[self.user.id, ],
                countdown=settings.POPOVERS_DISPLAY_REFRESH_PERIOD)
            self.task_id = res.id
            self.save()

    def update_display(self):
        upos = list(UserPopOver.objects
                    .filter(user=self.user, status=NEW))
        num_anchored_added = 0
        num_modal_added = 0
        with transaction.atomic():
            for upo in upos:
                if settings.POPOVERS_DISPLAY_HIGHEST_PRIORITY \
                        and upo.popover.priority == 1:
                    upo.display = DISPLAY
                elif upo.popover.type == MODAL \
                        and num_modal_added < settings.POPOVERS_DISPLAY_NEW_MODAL:
                    upo.display = DISPLAY
                    num_modal_added += 1
                elif upo.popover.type == ANCHORED \
                        and num_anchored_added < settings.POPOVERS_DISPLAY_NEW_ANCHORED:
                    upo.display = DISPLAY
                    num_anchored_added += 1
                upo.save()

        # clear task_id
        self.task_id = ''
        self.save()
