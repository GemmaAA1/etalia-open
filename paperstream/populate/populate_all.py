"""Init database
"""

from subprocess import call

call("./manage.py populate publisher")
call("./manage.py populate journal all")
call("./manage.py populate consumer pubmed --name pubmed1")
call("./manage.py populate consumer elsevier --name elsevier1")
call("./manage.py populate consumer arxiv --name elsevier1")