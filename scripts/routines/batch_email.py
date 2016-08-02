# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.conf import settings
from django.template.loader import get_template
from django.core.mail import EmailMultiAlternatives
from django.template import Context
from time import sleep

from_email = 'nicolas.pannetier@gmail.com'
email_file = 'beta_emails.txt'


def send_invite(to):
    to = [to, ]
    subject = 'An invitation to try Etalia'
    ctx = {'bucket_url': settings.EMAIL_STATIC_BUCKET,
           'root_url': 'etalia.io'}
    text_content = ''
    html_content = get_template(settings.INVITE_EMAIL_TEMPLATE)\
        .render(Context(ctx))
    email = EmailMultiAlternatives(subject,
                                   text_content,
                                   to=to,
                                   from_email=from_email)
    email.attach_alternative(html_content, "text/html")
    # email.content_subtype = 'html'
    email.send()


if __name__ == '__main__':
    with open(email_file, 'r') as f:
        emails = f.readlines()

    emails = [email.strip() for email in emails]

    for email in emails:
        sleep(0.5)
        send_invite(to=email)
