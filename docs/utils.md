### Show list of urls

    ./manage.py show_urls


### Create admin user

    ./manage.py shell
    from etalia.users.models import User
    User.objects.create_superuser(email='username@etalia.io', password='password')

### Shell plus (enhanced shell from django-extension)

    ./manage.py shell_plus 