from django.db import models

# Create your models here.

from django.db import models
from paperstream.core.models import TimeStampedModel
# Create your models here.


class EmailModel(TimeStampedModel):

    email = models.EmailField()

    def __str__(self):
        return '{email}'.format(email=self.email)
