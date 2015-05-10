from allauth.socialaccount.providers.oauth2.urls import default_urlpatterns
from .provider import MendeleyProvider

urlpatterns = default_urlpatterns(MendeleyProvider)

