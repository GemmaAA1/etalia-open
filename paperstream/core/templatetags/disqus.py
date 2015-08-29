# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django import template
import base64
import hashlib
import hmac
import json
import time

DISQUS_SECRET_KEY = 'eMWsm6qeNkDHzdvLViScWPldyDVnmvAz4U79YjsCelOu58XnRPelrUimqTrhGrRw'
DISQUS_PUBLIC_KEY = 'w2W0iBEJwGE49PjupwQxDnfzC9ayliEvctiGwbmVb63uHIXNTZLgreJDNRvvBOap'

register = template.Library()

@register.simple_tag
def get_disqus_sso(request):
    # create a JSON packet of our data attributes
    data = json.dumps({
        'id': request.user.id,
        'username': request.user.username,
        'email': request.user.email,
    })
    # encode the data to base64
    message = base64.b64encode(data.encode('utf-8')).decode()
    # generate a timestamp for signing the message
    timestamp = int(time.time())

    key = DISQUS_SECRET_KEY.encode('utf-8')
    msg = ('%s %s' % (message, timestamp)).encode('utf-8')
    digestmod = hashlib.sha1

    # generate our hmac signature
    sig = hmac.HMAC(key, msg, digestmod).hexdigest()


# return a script tag to insert the sso message
    return """<script type="text/javascript">
    var disqus_config = function() {
        this.page.remote_auth_s3 = "%(message)s %(sig)s %(timestamp)s";
        this.page.api_key = "%(pub_key)s";
    }
    </script>""" % dict(
        message=message,
        timestamp=timestamp,
        sig=sig,
        pub_key=DISQUS_PUBLIC_KEY,
    )