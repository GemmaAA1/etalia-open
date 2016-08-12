# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.conf import settings
from django.template.loader import get_template
from django.core.mail import EmailMultiAlternatives
from django.template import Context
from .models import UserInvited


def send_invite_email(email_to=None,
                      root_url=None,
                      on_behalf=None):
    subject = 'An invitation to try Etalia'
    to = [email_to]
    from_email = 'nicolas.pannetier@etalia.io'
    if on_behalf:
        user_name= '{first} {last}'.format(
                first=on_behalf.first_name,
                last=on_behalf.last_name
        )
    else:
        user_name = None
    ctx = {
        'bucket_url': settings.EMAIL_STATIC_BUCKET,
        'root_url': root_url,
        'on_behalf': user_name,
    }
    text_content = ''
    html_content = get_template(settings.INVITE_EMAIL_TEMPLATE)\
        .render(Context(ctx))
    email = EmailMultiAlternatives(subject, text_content, to=to, from_email=from_email)
    email.attach_alternative(html_content, "text/html")
    # email.content_subtype = 'html'
    email.send()

    # save to database
    UserInvited.objects.create(from_user=on_behalf, email_to=email_to)
