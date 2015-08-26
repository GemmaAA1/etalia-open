from .base import *
from core.utils import get_env_variable

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = get_env_variable('PAP_SECRET_KEY')