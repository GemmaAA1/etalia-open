from django import forms

from .models import Affiliation


class SocialPrimaryForm(forms.Form):

    first_name = forms.RegexField(
        regex=r'^[a-zA-Z\s]+$',
        max_length=255,
        label="First Name",
        error_messages={'invalid': "First name must contain only letters"},
        widget=forms.TextInput(attrs={'placeholder': 'First Name'})
    )

    last_name = forms.RegexField(
        regex=r'^[a-zA-Z\s]+$',
        max_length=255,
        label="Last Name",
        error_messages={'invalid': "Last name must contain only letters"},
        widget=forms.TextInput(attrs={'placeholder': 'Last Name'})
    )

    email = forms.EmailField(
        max_length=255,
        label="Email",
        widget=forms.TextInput(attrs={'placeholder': 'Email'})
    )


class CustomSignupForm(forms.Form):

    first_name = forms.RegexField(
        regex=r'^[a-zA-Z\s]+$',
        max_length=30,
        label="First Name",
        error_messages={'invalid': "First name must contain only letters"},
        widget=forms.TextInput(attrs={'placeholder': 'First Name'})
    )

    last_name = forms.RegexField(
        regex=r'^[a-zA-Z\s]+$',
        max_length=30,
        label="Last Name",
        error_messages={'invalid': "Last name must contain only letters"},
        widget=forms.TextInput(attrs={'placeholder': 'Last Name'})
    )

    def signup(self, request, user):
        user.first_name = self.cleaned_data['first_name']
        user.last_name = self.cleaned_data['last_name']
        user.email = self.cleaned_data['email']
        user.save()

        # Asynchro Init library user
        # TODO: Is there a way to start this earlier ?
        if hasattr(self, 'sociallogin'):
            pass


class AffiliationForm(forms.ModelForm):

    class Meta:
        model = Affiliation
        fields = ('department', 'institution', 'city', 'state', 'country')