from django.db import models

# Create your models here.

from django.db import models

# Create your models here.


class TimeStampedModel(models.Model):
    """
    An abstract base class model that provides selfupdating
    ``created`` and ``modified`` fields.
    """
    created = models.DateTimeField(auto_now_add=True)
    modified = models.DateTimeField(auto_now=True)

    class Meta:
        abstract = True


class EmailModel(TimeStampedModel):

    email = models.EmailField()

    def __str__(self):
        return '{email}'.format(email=self.email)
