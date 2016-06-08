### After rebase/merge :

* Check diff between common.py and common.py.dist
* $ ./manage.py migrate
* $ ./scripts/routines/update.py

### Show list of urls

$ ./manage.py show_urls


### Create admin user
./manage.py shell
from etalia.users.models import User
User.objects.create_superuser(email='username@etalia.io', password='password')