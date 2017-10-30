from django.contrib import admin

# Register your models here.

from django.contrib import admin

from .models import PopOver, UserPopOver

admin.site.register(PopOver)
admin.site.register(UserPopOver)
