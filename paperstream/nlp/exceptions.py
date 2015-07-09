

class ModelError(Exception):
    """Top-level error for model exception.

    This exception should normally not be raised, only subclasses of this
    exception."""

    def __str__(self):
        return getattr(self, 'message', '')


class StatusError(ModelError):
    """The state of the model is wrong

    This generally means that the model is currently doing something else
    and/or the current task cannot be processed.
    """

    message = 'You can do this in the current model status'