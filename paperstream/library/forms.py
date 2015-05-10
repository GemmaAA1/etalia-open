from django.forms import ValidationError, ModelForm

from .models import Journal, Paper, Author


class JournalForm(ModelForm):

    class Meta:
        model = Journal


class AuthorForm(ModelForm):

    class Meta:
        model = Author


class PaperForm(ModelForm):

    class Meta:
        model = Paper
