from django.shortcuts import render

# Create your views here.

def view_library(request):
    render(request, 'library.html')

def view_papers(request, paper_id):
    pass

def view_journals(request, journal_id):
    pass
