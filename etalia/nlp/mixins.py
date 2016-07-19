# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import boto
import glob
import tarfile
import logging
import lz4
from boto.s3.key import Key
from progressbar import ProgressBar, Percentage, Bar, ETA


class S3Mixin(object):

    BUCKET_NAME = ''
    AWS_ACCESS_KEY_ID = ''
    AWS_SECRET_ACCESS_KEY = ''
    PATH = ''
    name = ''

    def s3callback(self, complete, tot):
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

            # Compressing
            logging.info('{} LZ4 compressing...'.format(self.name))
            tar_path = os.path.join(self.PATH, '{}.tar'.format(self.name))
            lz4_path = os.path.join(self.PATH, '{}.tar.lz4'.format(self.name))

            with tarfile.open(tar_path, 'w') as in_file:
                for filename in glob.glob(os.path.join(self.PATH,
                                                       '{stem}.*'.format(
                                                           stem=self.name))):
                    in_file.add(filename, arcname=os.path.split(filename)[1])

            with open(tar_path, 'rb') as in_file:
                with open(lz4_path, 'wb') as out_file:
                    out_file.write(lz4.compress(in_file.read()))

            # Upload
            logging.info('Uploading {} on s3...'.format(self.name))
            # create a key to keep track of our file in the storage
            k = Key(bucket)
            k.key = '{}.tar.lz4'.format(self.name)
            k.set_contents_from_filename(lz4_path, cb=self.s3callback, num_cb=100)
            # remove tar and lz4 file
            os.remove(tar_path)
            os.remove(lz4_path)
            logging.info('Upload {} DONE'.format(self.name))

        except Exception:
            raise

    def pull_from_s3(self):
        """Download object from Amazon s3 bucket"""
        try:
            bucket_name = self.BUCKET_NAME
            conn = boto.connect_s3(self.AWS_ACCESS_KEY_ID,
                                   self.AWS_SECRET_ACCESS_KEY)
            bucket = conn.get_bucket(bucket_name)
            key = '{}.tar.lz4'.format(self.name)
            item = bucket.get_key(key)

            # Downloading
            lz4_path = os.path.join(self.PATH, key)
            logging.info('Downloading {0} from s3 to {1}...'.format(self.name,
                                                              self.PATH))
            item.get_contents_to_filename(lz4_path,
                                          cb=self.s3callback,
                                          num_cb=100)

            # Decompressing
            logging.info('{} Decompressing...'.format(self.name))
            tar_path = os.path.join(self.PATH, '{}.tar'.format(self.name))
            with open(lz4_path, 'rb') as in_file:
                with open(tar_path, 'wb') as out_file:
                    out_file.write(lz4.decompress(in_file.read()))

            with tarfile.open(tar_path, 'r') as in_file:
                    in_file.extractall(self.PATH)
            # remove tar file
            try:
                os.remove(lz4_path)
            except FileNotFoundError:
                pass
            try:
                os.remove(tar_path)
            except FileNotFoundError:
                pass
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
            key = self.name + '.tar.lz4'
            item = bucket.get_key(key)
            if item:
                logging.info('{} Deleting on s3...'.format(self.name))
                item.delete()
                logging.info('{} Done'.format(self.name))
        except Exception:
            raise
