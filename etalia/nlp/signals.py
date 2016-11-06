# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.db.models.signals import post_delete
from django.dispatch import receiver

from etalia.library.models import Paper
from etalia.threads.models import Thread


@receiver(post_delete, sender=Paper)
def remove_entry_in_paperengine(sender, instance, using, **kwargs):
    """Remove entry from PaperEngine data structure"""
    from .tasks import pe_dispatcher
    # Broadcast task to all PaperEngine
    pe_dispatcher.apply_async(args=['remove_entry_id', instance.id],
                              queue='broadcast_pe_tasks')


@receiver(post_delete, sender=Thread)
def remove_entry_in_threadengine(sender, instance, using, **kwargs):
    """Remove entry from ThreadEngine data structure"""
    from .tasks import te_dispatcher
    # Broadcast task to all ThreadEngine
    te_dispatcher.apply_async(args=['remove_entry_id', instance.id],
                              queue='broadcast_te_tasks')
