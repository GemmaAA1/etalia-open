from django.shortcuts import render, redirect
from django.core.urlresolvers import reverse

def home(request):
    if request.user.is_authenticated():
        return redirect('feeds:feed')
    else:
        return render(request, 'landing.html')

def test(request):
    return render(request, 'test.html', {})


class TitleSearchMixin(object):

    def get_queryset(self):
        # Fetch the queryset from the parent's get_queryset
        queryset = super(TitleSearchMixin, self).get_queryset()
        # Get the q GET parameter
        q = self.request.GET.get("q")
        if q:
            # return a filtered queryset
            return queryset.filter(title__icontains=q)
        # No q is specified so we return queryset
        return queryset