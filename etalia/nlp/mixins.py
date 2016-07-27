# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import os
import boto3
import botocore
import glob
import tarfile
import logging
import lz4
import threading
from progressbar import ProgressBar, Percentage, Bar, ETA
from django.core.exceptions import ImproperlyConfigured


class S3Mixin(object):

    BUCKET_NAME = ''
    AWS_ACCESS_KEY_ID = ''
    AWS_SECRET_ACCESS_KEY = ''
    PATH = ''
    name = ''

    class ProgressPercentage(object):

        def __init__(self, content_length=None):
            self._size = content_length
            self._seen_so_far = 0
            self._pbar = None

        def __call__(self, bytes_amount):
            self._seen_so_far += bytes_amount
            percentage = (self._seen_so_far / self._size) * 100
            if self._pbar is None:
                self._pbar = ProgressBar(widgets=[Percentage(), Bar(), ' ', ETA()],
                                         maxval=100, term_width=75).start()
            if not self._pbar.finished:
                self._pbar.update(percentage)
                if percentage == 100:
                    # close progress bar
                    self._pbar.finish()

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

    def connect_and_check_bucket(self):

        bucket_name = self.BUCKET_NAME
        s3 = boto3.resource(
            's3',
            aws_access_key_id=self.AWS_ACCESS_KEY_ID,
            aws_secret_access_key=self.AWS_SECRET_ACCESS_KEY)

        # Check existence of bucket
        try:
            s3.meta.client.head_bucket(Bucket=bucket_name)
        except botocore.exceptions.ClientError as e:
            # If a client error is thrown, then check that it was a 404 error.
            # If it was a 404 error, then the bucket does not exist.
            error_code = int(e.response['Error']['Code'])
            if error_code == 404 or 403:
                raise ImproperlyConfigured('bucket name ({}) does not exist'.format(bucket_name))

        return s3

    def get_key(self, compress=True):
        # Build key
        if compress:
            return '{}.tar.lz4'.format(self.name)
        else:
            return '{}.tar'.format(self.name)

    def get_keys(self):
        return [self.get_key(compress=False), self.get_key()]

    def push_to_s3(self, compress=True):
        """Upload object to Amazon s3 bucket"""

        s3 = self.connect_and_check_bucket()

        # Archive
        logging.info('Archiving {}...'.format(self.name))
        tar_path = os.path.join(self.PATH, '{}.tar'.format(self.name))
        lz4_path = os.path.join(self.PATH, '{}.tar.lz4'.format(self.name))

        with tarfile.open(tar_path, 'w') as in_file:
            for filename in glob.glob(os.path.join(self.PATH,
                                                   '{stem}.*'.format(
                                                       stem=self.name))):
                in_file.add(filename, arcname=os.path.split(filename)[1])

        # Compress
        if compress:
            logging.info('{} LZ4 compressing...'.format(self.name))
            with open(tar_path, 'rb') as in_file:
                with open(lz4_path, 'wb') as out_file:
                    out_file.write(lz4.compress(in_file.read()))
            upload_path = lz4_path
        else:
            upload_path = tar_path

        # Get key
        key = self.get_key(compress=compress)

        # Upload
        logging.info('Uploading {} on s3...'.format(self.name))
        s3.meta.client.upload_file(
            upload_path,
            self.BUCKET_NAME,
            key,
            Callback=self.ProgressPercentage(float(os.path.getsize(upload_path))))

        # Clean
        os.remove(tar_path)
        if compress:
            os.remove(lz4_path)
        logging.info('Upload {} DONE'.format(self.name))

    def pull_from_s3(self, compress=True):
        """Download object from Amazon s3 bucket"""

        s3 = self.connect_and_check_bucket()

        # Get key
        key = self.get_key(compress=compress)

        # Download
        download_path = os.path.join(self.PATH, key)
        logging.info('Downloading {0} from s3 to {1}...'.format(self.name,
                                                          self.PATH))
        content_length = s3.Object(self.BUCKET_NAME, key).content_length
        s3.meta.client.download_file(
            self.BUCKET_NAME,
            key,
            download_path,
            Callback=self.ProgressPercentage(content_length))

        # Decompress
        logging.info('{} Decompressing...'.format(self.name))
        tar_path = os.path.join(self.PATH, '{}.tar'.format(self.name))
        if compress:
            with open(download_path, 'rb') as in_file:
                with open(tar_path, 'wb') as out_file:
                    out_file.write(lz4.decompress(in_file.read()))

        # Expand
        with tarfile.open(tar_path, 'r') as in_file:
                in_file.extractall(self.PATH)

        # Clean
        try:
            os.remove(download_path)
        except FileNotFoundError:
            pass
        try:
            os.remove(tar_path)
        except FileNotFoundError:
            pass
        logging.info('{} Done'.format(self.name))

    def delete_on_s3(self):
        """Delete object from Amazon s3 bucket"""

        s3 = self.connect_and_check_bucket()
        bucket = s3.Bucket(self.BUCKET_NAME)

        # Get keys
        keys = self.get_keys()

        # Delete
        for obj in bucket.objects.all():
            if obj.key in keys:
                obj.delete()
