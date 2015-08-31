# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

from django.contrib import admin

# Register your models here.


from .models import Paper, Author, Journal, Publisher



class AuthorAdmin(admin.ModelAdmin):
    search_fields = ('title', )
    list_display = ('email', 'first_name', 'last_name')

admin.site.register(Author, AuthorAdmin)


class PaperAdmin(admin.ModelAdmin):
    search_fields = ('title', 'author', )
    list_display = ('short_title', 'print_compact_first_author', 'journal',
                    'date_ep', 'date_pp')

admin.site.register(Paper, PaperAdmin)


class JournalAdmin(admin.ModelAdmin):
    search_fields = ('title', )
    list_display = ('title', 'publisher', 'id_issn', 'id_eissn', 'id_arx',
                    'id_oth', 'period', 'url')

admin.site.register(Journal, JournalAdmin)

admin.site.register(Publisher)