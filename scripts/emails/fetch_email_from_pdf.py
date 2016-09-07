# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import

import requests
import PyPDF2
import re
import os

DOWNLOAD_DIR = '/tmp/pdfs'
LOG_FILE = 'log.txt'


def download_file(url):
    local_filename = os.path.join(DOWNLOAD_DIR, url.split('/')[-1])
    dir_name = os.path.dirname(local_filename)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
    return local_filename


def retrieve(file_name):
    # Define variables

    file_name_emails = '{base}_emails.txt'.format(base=os.path.splitext(file_name)[0])

    if not os.path.exists(LOG_FILE):
        with open(LOG_FILE, 'w+') as file:
            pass
        log_urls = []
    else:
        with open(LOG_FILE, 'rb') as file:
            log_urls = file.readlines()

    # Open urls.txt
    with open(file_name) as file:
        for url in file:
            url = url.rstrip()

            if url not in log_urls:

                print('processing {0}'.format(url))
                emails_found = []
                try:
                    # Download pdf
                    pdf_file_name = download_file(url)

                    # Read pdf
                    with open(pdf_file_name, 'rb') as pdf_file:
                        pdf_reader = PyPDF2.PdfFileReader(pdf_file)
                        for p in range(min(pdf_reader.numPages, 3)):
                            page_obj = pdf_reader.getPage(p)
                            textfile = page_obj.extractText()

                            # Regular expression to find email pattern
                            mailsrch = re.compile(r'[a-zA-Z0-9+_\-\.]+@[0-9a-zA-Z]+[\.-0-9a-zA-Z]*\.[a-zA-Z]+')
                            # Search and extract all email adresses from the PDD text files
                            res = mailsrch.findall(textfile)
                            if res:
                                for i in res:
                                    emails_found.append(i.lower())

                            # Try {xxx, yyy}@domain.zzz format sometimes found in latex
                            mailsrch2 = re.compile(r'\nf\n[\,a-zA-Z0-9+_\-\.\s]+\ng\n@[0-9a-zA-Z]+[\.-0-9a-zA-Z]*\.[a-zA-Z]+')
                            res = mailsrch2.findall(textfile)
                            if res:
                                for i in res:
                                    tmp = i.split('@')
                                    domain = tmp[1]
                                    names = tmp[0]
                                    names = names.rstrip('\ng\n')
                                    names = names.lstrip('\nf\n')
                                    names = names.split(',')
                                    names = [name.rstrip().lstrip() for name in names]
                                    for name in names:
                                        emails_found.append('{0}@{1}'.format(name, domain).lower())

                    # remove pdf
                    os.remove(pdf_file_name)
                except Exception:
                    pass

                # Write emails to file "email_exports/emails.txt"
                if emails_found:
                    with open(file_name_emails, 'w+') as file:
                        for email in emails_found:
                            file.write('{0}\n'.format(email))

                # Log url
                with open(LOG_FILE, 'w+') as file:
                    file.write('{0}\n'.format(url))
