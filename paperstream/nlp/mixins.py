# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from progressbar import ProgressBar, Percentage, Bar, ETA


class S3ProgressBarMixin(object):

    def callback(self, complete, tot):
        if not hasattr(self, 'pbar'):
            setattr(self, 'pbar', None)

        if not self.pbar:
            self.pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                                    maxval=tot).start()
        self.pbar.update(complete)
        if complete == tot:
            # close progress bar
            self.pbar.finish()
            self.pbar = None
