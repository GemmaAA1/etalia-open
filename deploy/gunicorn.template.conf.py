import multiprocessing

name = "{{ SITENAME }}"
bind = "unix:/tmp/{{ SITENAME }}.socket"
workers = multiprocessing.cpu_count() * 2 + 1
user = "{{ USER }}"
