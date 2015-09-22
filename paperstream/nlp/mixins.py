# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import boto
import glob
import tarfile
import logging
from boto.s3.key import Key
from progressbar import ProgressBar, Percentage, Bar, ETA


class S3Mixin(object):

    BUCKET_NAME = ''
    AWS_ACCESS_KEY_ID = ''
    AWS_SECRET_ACCESS_KEY = ''
    PATH = ''

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

    def push_to_s3(self):
        """Upload object to Amazon s3 bucket"""
        try:
            bucket_name = self.BUCKET_NAME
            conn = boto.connect_s3(self.AWS_ACCESS_KEY_ID,
                                   self.AWS_SECRET_ACCESS_KEY)
            bucket = conn.get_bucket(bucket_name)

            # Compress file
            key = '{}.tar.gz'.format(self.name)
            tar_name = os.path.join(self.PATH, key)
            tar = tarfile.open(tar_name, 'w:gz')
            logging.info('{} Compressing...'.format(self.name))
            for filename in glob.glob(os.path.join(self.PATH,
                                                   '{0}.mod*'.format(self.name))):
                tar.add(filename, arcname=os.path.split(filename)[1])
            tar.close()
            logging.info('↑ {} Uploading...'.format(self.name))

            # create a key to keep track of our file in the storage
            k = Key(bucket)
            k.key = key
            k.set_contents_from_filename(tar_name, cb=self.callback, num_cb=100)
            # remove tar file
            os.remove(tar_name)
            logging.info('↑ {} Upload DONE'.format(self.name))

        except Exception:
            raise

    def download_from_s3(self):
        """Download object from Amazon s3 bucket"""
        try:
            bucket_name = self.BUCKET_NAME
            conn = boto.connect_s3(self.AWS_ACCESS_KEY_ID,
                                   self.AWS_SECRET_ACCESS_KEY)
            bucket = conn.get_bucket(bucket_name)
            key = self.name + '.tar.gz'
            item = bucket.get_key(key)
            tar_path = os.path.join(self.PATH, key)
            logging.info('↓ {} Downloading...'.format(self.name))
            item.get_contents_to_filename(tar_path,
                                          cb=self.callback,
                                          num_cb=100)
            logging.info('{} Decompressing...'.format(self.name))
            tar = tarfile.open(tar_path, 'r:gz')
            tar.extractall(self.PATH)
            tar.close()
            # remove tar file
            os.remove(tar_path)
            logging.info('{} Done'.format(self.name))
        except Exception:
            raise

    def delete_on_s3(self):
        """Delete object from Amazon s3 bucket"""
        try:
            bucket_name = self.BUCKET_NAME
            conn = boto.connect_s3(self.AWS_ACCESS_KEY_ID,
                                   self.AWS_SECRET_ACCESS_KEY)
            bucket = conn.get_bucket(bucket_name)
            key = self.name + '.tar.gz'
            item = bucket.get_key(key)
            if item:
                logging.info('{} Deleting...'.format(self.name))
                item.delete()
                logging.info('{} Done'.format(self.name))
        except Exception:
            raise
