# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
from anymail.message import attach_inline_image_file
from django.conf import settings
from django.template.loader import render_to_string
from django.core.mail import EmailMultiAlternatives


class Email(object):
    """Class that build and send email based on templates"""

    IMG_DIR = str(settings.ROOT_DIR.path('etalia/templates/emails/img/'))
    ROOT_URL = 'https://etalia.io'

    def __init__(self, *args, **kwargs):
        self.email = None
        self.template = kwargs.get('template', '')
        self.cids = kwargs.get('cids', {})
        self.tags = kwargs.get('tags', [])
        self.metadata = kwargs.get('metadata', {})
        self.subject = kwargs.get('subject', '')
        self.from_email = kwargs.get('from_email', '')
        self.reply_to = kwargs.get('reply_to', [])
        self.to = kwargs.get('to', [])
        self.cc = kwargs.get('cc', [])
        self.bcc = kwargs.get('bcc', [])
        self.extra_ctx = kwargs.get('extra_ctw', {})
        self.img_dir = kwargs.get('img_dir', self.IMG_DIR)
        self.root_url = kwargs.get('root_url', self.ROOT_URL)
        self.campaign_id = kwargs.get('campaign_id')
        self.build()

    def get_context(self):
        ctx = {'root_url': self.root_url}
        ctx.update(self.extra_ctx)
        return ctx

    def update_context_with_images(self, ctx):
        ctx.update(self.get_inline_images())
        return ctx

    def get_plaintext_template(self):
        return self.template.split('-')[0] + '-plain.txt'

    def get_inline_images(self):
        """Process inline images"""
        inline_images = {}
        for cid, path in self.cids.items():
            inline_images[cid] = attach_inline_image_file(
                self.email,
                os.path.join(self.img_dir, path),
            )
        return inline_images

    def get_metadata(self, **kwargs):
        return self.metadata

    def build(self):
        ctx = self.get_context()
        body = render_to_string(self.get_plaintext_template(), ctx)
        self.email = EmailMultiAlternatives(subject=self.subject,
                                            body=body,
                                            from_email=self.from_email,
                                            reply_to=self.reply_to,
                                            to=self.to,
                                            cc=self.cc,
                                            bcc=self.bcc)
        ctx = self.update_context_with_images(ctx)
        html_content = render_to_string(self.template, ctx)
        self.email.attach_alternative(html_content, "text/html")
        self.email.tags = self.tags
        self.email.campaign = self.campaign_id
        self.email.metadata = self.get_metadata()
        self.email.track_clicks = True

        return self.email

    def send_to(self, to):
        self.email.to = to
        self.send()

    def send(self):
        self.email.send()
