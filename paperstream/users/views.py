from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect
from django.db.models import Q
from django.contrib.auth import logout as auth_logout, login
from django.views.generic import UpdateView, FormView
from django.views.generic.list import ListView
from django.conf import settings
from django.contrib.auth import get_user_model
from django.http import JsonResponse

from braces.views import LoginRequiredMixin

from paperstream.library.models import Paper
from paperstream.core.mixins import AjaxableResponseMixin

from .forms import UserBasicForm, UserAffiliationForm, UpdateUserBasicForm, \
    UserAuthenticationForm, UserSettingsForm
from .models import Affiliation
from .tasks import update_lib as async_update_lib

User = get_user_model()

def logout(request):
    """Logs out user"""
    auth_logout(request)
    return redirect('/')

def done(request):
    redirect('landing')

# User authentication
# ---------------
class UserLoginView(AjaxableResponseMixin, FormView):

    form_class = UserAuthenticationForm
    redirect_field_name = settings.LOGIN_REDIRECT_URL

    def form_valid(self, form):
        login(self.request, form.get_user())
        response = super(AjaxableResponseMixin, self).form_valid(form)
        if self.request.is_ajax():
            data = self.get_ajax_data()
            return JsonResponse(data)
        else:
            return response

    def get_success_url(self):
        return reverse('feeds:main')

    def get_ajax_data(self):
        data = {'email': self.request.user.email,
                'redirect': self.get_success_url(),
                }
        return data

ajax_signin = UserLoginView.as_view()


class UserBasicInfoSignupView(AjaxableResponseMixin, FormView):

    form_class = UserBasicForm
    template_name = 'user/basic_info.html'

    def get_success_url(self, **kwargs):
        return reverse('social:complete',
        kwargs={'backend': self.request.session['partial_pipeline']['backend']})

    def get_context_data(self, **kwargs):
        context = super(UserBasicInfoSignupView, self).get_context_data(**kwargs)
        context['stage'] = 'basic_info'
        return context

    def get_initial(self):
        initial = super(UserBasicInfoSignupView, self).get_initial()
        details = self.request.session['partial_pipeline']['kwargs']['details']
        initial['first_name'] = details.get('first_name', '')
        initial['last_name'] = details.get('last_name', '')
        initial['email'] = details.get('email', '')
        return initial

    def form_valid(self, form):
        self.request.session['basic_info'] = {
            'first_name': form.cleaned_data['first_name'],
            'last_name': form.cleaned_data['last_name'],
            'email': form.cleaned_data['email'],
            }
        return super(UserBasicInfoSignupView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'redirect': self.get_success_url()}
        return data

require_basic_info = UserBasicInfoSignupView.as_view()


class UserAffiliationSignupView(FormView):

    form_class = UserAffiliationForm
    template_name = 'user/basic_info.html'

    def get_success_url(self, **kwargs):
        return reverse('social:complete',
        kwargs={'backend': self.request.session['partial_pipeline']['backend']})

    def get_context_data(self, **kwargs):
        context = super(UserAffiliationSignupView, self).get_context_data(**kwargs)
        context['stage'] = 'affiliation'
        return context

    def get_initial(self):
        initial = super(UserAffiliationSignupView, self).get_initial()
        details = self.request.session['partial_pipeline']['kwargs']['details']
        return dict(details.get('tmp_affiliation', ''), **initial)

    def form_invalid(self, form):
        if Affiliation.objects.filter(**form.cleaned_data).exists():  # invalid because of uniqueness
            affiliation = Affiliation.objects.get(**form.cleaned_data)
            self.request.session['affiliation_pk'] = affiliation.pk
            return super(UserAffiliationSignupView, self).form_valid(form)
        else:
            return super(UserAffiliationSignupView, self).form_invalid(form)

    def form_valid(self, form):
        affiliation = form.save()
        self.request.session['affiliation_pk'] = affiliation.pk
        return super(UserAffiliationSignupView, self).form_valid(form)

require_affiliation = UserAffiliationSignupView.as_view()


def validation_sent(request):
    context = {
        'stage': 'validation_sent',
        'email': request.session.get('email_validation_address')
    }
    return render(request, 'user/basic_info.html', context)


# User Library
# ---------------
class UserLibraryView(LoginRequiredMixin, ListView):
    model = Paper
    template_name = 'user/library.html'
    paginate_by = settings.ITEMS_PER_PAGE
    context_object_name = 'paper_list'

    def get_context_data(self, **kwargs):
        context = super(UserLibraryView, self).get_context_data(**kwargs)
        return context

    def get_queryset(self):
        query_set = self.request.user.lib.papers.all()
        return self.filter_queryset(query_set)

    def filter_queryset(self, queryset):
        # Get the q GET parameter
        q = self.request.GET.get("query")
        if q:
            # return a filtered queryset
            return queryset.filter(Q(title__icontains=q) |
                                   Q(abstract__icontains=q) |
                                   Q(journal__title__icontains=q) |
                                   Q(authors__last_name__icontains=q) |
                                   Q(authors__first_name__icontains=q)).distinct()
        # No q is specified so we return queryset
        return queryset

library = UserLibraryView.as_view()


# User profile update
# -------------------
class UpdateUserBasicInfoView(LoginRequiredMixin, AjaxableResponseMixin,
                              UpdateView):
    form_class = UpdateUserBasicForm

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        return reverse('core:home')

    def get_ajax_data(self):
        data = {'first_name': self.request.user.first_name,
                'last_name': self.request.user.last_name}
        return data

ajax_update_basic_info = UpdateUserBasicInfoView.as_view()


class UserAffiliationUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                                FormView):
    form_class = UserAffiliationForm

    def get_object(self, queryset=None):
        return self.request.user.affiliation

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        if Affiliation.objects.filter(**form.cleaned_data).exists():
            self.request.user.affiliation = \
                Affiliation.objects.get(**form.cleaned_data)
            self.request.user.save()
            return super(UserAffiliationUpdateView, self).form_valid(form)
        else:
            return super(UserAffiliationUpdateView, self).form_invalid(form)

    def form_valid(self, form):
        affiliation = form.save()
        self.request.user.affiliation = affiliation
        self.request.user.save()
        return super(UserAffiliationUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'department': self.request.user.affiliation.department,
                'institution': self.request.user.affiliation.institution,
                'city': self.request.user.affiliation.city,
                'state': self.request.user.affiliation.state,
                'country': self.request.user.affiliation.country}
        return data

ajax_update_affiliation = UserAffiliationUpdateView.as_view()


class UserSettingsUpdateView(LoginRequiredMixin, AjaxableResponseMixin,
                             FormView):
    form_class = UserSettingsForm

    def get_object(self, queryset=None):
        return self.request.user.settings

    def get_success_url(self):
        return reverse('core:home')

    def form_invalid(self, form):
        return super(UserSettingsUpdateView, self).form_invalid(form)

    def form_valid(self, form):
        self.request.user.settings.time_lapse = form.cleaned_data['time_lapse']
        self.request.user.settings.model = form.cleaned_data['model']
        self.request.user.settings.scoring_method = \
            form.cleaned_data['scoring_method']
        self.request.user.settings.save()
        return super(UserSettingsUpdateView, self).form_valid(form)

    def get_ajax_data(self):
        data = {'model': self.request.user.settings.model,
                'time_lapse': self.request.user.settings.time_lapse,
                'scoring_method': self.request.user.settings.scoring_method,
                }
        return data

ajax_update_settings = UserSettingsUpdateView.as_view()

@login_required
def ajax_user_lib_count_papers(request):
    if request.method == 'GET':
        if request.user.lib.state == 'IDL':
            data = {'done': True}
        else:
            data = {'done': False,
                    'message': request.user.lib.count_papers}
        return JsonResponse(data)

@login_required
def async_update_user_lib(request):
    user = request.user
    provider_name = user.social_auth.first().provider
    async_update_lib.apply_async(args=[user.pk, provider_name])
    request.user.lib.set_state('ING')
    return redirect('feeds:main')

