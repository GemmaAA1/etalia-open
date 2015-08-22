import hashlib
from django.conf import settings

def get_disqus_username(user_id):
    m = hashlib.md5(user_id)
    hashed_id = m.hexdigest()
    return '{short_name}-{hashed_id}'.format(
        short_name=settings.DISQUS_WEBSITE_SHORTNAME,
        hashed_id=hashed_id)
