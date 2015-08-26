from .base import *
from core.utils import get_env_variable

SECRET_KEY = get_env_variable('PAP_SECRET_KEY')