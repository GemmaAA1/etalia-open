Easy hacks to tweak data from django shell
========================

```bash
./manage.py shell
```

## Get user

```python
from django.contrib.auth import get_user_model
User = get_user_model()
user = User.objects.get(email='trucfortest@gmail.com')
```

## Update user initialization steps (require user)


```python
user.init_step = 'STR'      # choices are 'TRE' for trend, 'LIB' for library, 'IDL' for idle
user.save()
```

## Update user stream state (require user)


```python
stream = user.streams.first()
stream.set_state('ING')     # choices are 'ING' for processING, or 'IDL' for idle
```

## Update user stream state (require user)


```python
trend = user.trends.first()
trend.set_state('ING')     # choices are 'ING' for processING, or 'IDL' for idle
```