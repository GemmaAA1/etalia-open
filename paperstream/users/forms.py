from django import forms
from django.contrib.auth import get_user_model
from django.contrib.auth.forms import AuthenticationForm
from .models import Affiliation, UserSettings
from .validators import validate_first_name, validate_last_name

from nlp.models import Model

User = get_user_model()


class UserAuthenticationForm(AuthenticationForm):

    error_messages = {
        'invalid_login': "Please enter a correct %(username)s and password.",
        'inactive': "This account is inactive.",
    }


class UserBasicForm(forms.ModelForm):

    error_messages = {
        'password_mismatch': "The two password fields didn't match.",
    }

    password1 = forms.CharField(label="Password",
        widget=forms.PasswordInput,
        initial='')

    password2 = forms.CharField(label="Password confirmation",
        widget=forms.PasswordInput,
        help_text="Enter the same password as above, for verification.",
        initial='')

    def __init__(self, *args, **kwargs):
        super(UserBasicForm, self).__init__(*args, **kwargs)
        self.fields['first_name'].validators.append(validate_first_name)
        self.fields['last_name'].validators.append(validate_last_name)

    def clean_first_name(self):
        first_name = self.cleaned_data['first_name']
        return first_name.strip()

    def clean_last_name(self):
        last_name = self.cleaned_data['last_name']
        return last_name.strip()

    def clean_password2(self):
        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError(
                self.error_messages['password_mismatch'],
                code='password_mismatch',
            )
        return password2

    class Meta:
        model = User
        fields = ('first_name', 'last_name', 'email')


class UpdateUserBasicForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UpdateUserBasicForm, self).__init__(*args, **kwargs)
        self.fields['first_name'].validators.append(validate_first_name)
        self.fields['last_name'].validators.append(validate_last_name)

    def clean_first_name(self):
        first_name = self.cleaned_data['first_name']
        return first_name.strip()

    def clean_last_name(self):
        last_name = self.cleaned_data['last_name']
        return last_name.strip()

    class Meta:
        model = User
        fields = ('first_name', 'last_name')


class UserAffiliationForm(forms.ModelForm):

    class Meta:
        model = Affiliation
        fields = ('department', 'institution', 'city', 'state', 'country')


class UserProfileForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ('first_name', 'last_name', 'email', 'affiliation')


class UserSettingsForm(forms.ModelForm):

    model = forms.ModelChoiceField(queryset=Model.objects.all(),
                                   to_field_name='name')

    def __init__(self, *args, **kwargs):
        super(UserSettingsForm, self).__init__(*args, **kwargs)
        # self.fields['model'].queryset = Model.objects.all()

    class Meta:
        model = UserSettings
        fields = ('model', 'time_lapse')

