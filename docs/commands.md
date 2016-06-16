### After rebase/merge :

* Check diff between common.py and common.py.dist
* $ ./manage.py migrate
* $ ./scripts/routines/update.py

### Show list of urls

$ ./manage.py show_urls

## Install

* ./scripts/routines/update.py -i : (init) va peupler la db avec des data et initialiser les objet de NLP
* ./scripts/routines/update.py -l : (load) essaye de loader les data partir de la fixture json


### Create admin user
./manage.py shell
from etalia.users.models import User
User.objects.create_superuser(email='username@etalia.io', password='password')
