# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import tarfile
import boto
from boto.s3.key import Key
import os
import glob


bucket_name = 'paperstream-nlp-models'
conn = boto.connect_s3('AKIAJP4QVWCJZTBCDW7A','np3BxaZhtxAp1i9pYQ6g1lEvb5KluUBR/DgisDu4')
bucket = conn.get_bucket(bucket_name)
NLP_MODELS_PATH = '/Users/nicolaspannetier/Projects/paperstream/paperstream_project/nlp_data/mods'

key = 'test.tar.gz'
tar_name = os.path.join(NLP_MODELS_PATH, key)
tar = tarfile.open(tar_name, 'w:gz')
for filename in glob.glob(os.path.join(NLP_MODELS_PATH, 'test*')):
    tar.add(filename, arcname=os.path.split(filename)[1])
tar.close()
# create a key to keep track of our file in the storage
k = Key(bucket)
k.key = key
k.set_contents_from_filename(tar_name)
# remove tar file
os.remove(tar_name)


item = bucket.get_key(key)
tar_path = os.path.join(NLP_MODELS_PATH, key)
item.get_contents_to_filename(tar_path)
tar = tarfile.open(tar_path, 'r:gz')
tar.extractall()
# remove tar file
os.remove(key)
