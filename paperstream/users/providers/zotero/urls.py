from allauth.socialaccount.providers.oauth.urls import default_urlpatterns

from .provider import ZoteroProvider

urlpatterns = default_urlpatterns(ZoteroProvider)
