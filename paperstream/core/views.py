from django.shortcuts import render, redirect


def home(request):
    # TODO: remove this and make url in template dynamic
    host = request.get_host()
    return render(request, 'landing.html', {'host': host})

def test(request):
    return render(request, 'test.html', {})