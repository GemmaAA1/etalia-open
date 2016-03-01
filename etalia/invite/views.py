from django.shortcuts import render, redirect
from django.http import JsonResponse
from django.core.mail import send_mail

# Create your views here.
from .models import EmailModel


def home(request):
    return render(request, 'invite/landing.html')


def request_invite(request):
    if request.POST:
        email = request.POST.get('email')
        # save email to db
        EmailModel.objects.create(email=email)

        # send email
        send_mail('Etalia - Invite Request', email, 'invite@etalia.io',
                  ['etalia.io@gmail.com'], fail_silently=False)

        return JsonResponse(data={'done': True})
    else:
        redirect('invite:home')
