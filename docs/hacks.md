Easy hacks to tweak data from django shell
========================

```bash
./manage.py shell_plus          # shell_plus is the enhanced shell from django-extension
```

## Get user

```python
user = User.objects.get(email='trucfortest@gmail.com')
```

## Update user initialization steps


```python
user.init_step = 'STR'      # choices are 'TRE' for TREnd, 'STR' for STReam, 'LIB' for LIBrary, 'IDL' for IDLe
user.save()
```

## Update user stream state 


```python
stream = user.streams.first()
stream.set_state('ING')     # choices are 'ING' for processING, or 'IDL' for IDLe
```

## Update user stream state


```python
trend = user.trends.first()
trend.set_state('ING')     # choices are 'ING' for processING, or 'IDL' for IDLe
```

## Delete all users


```python
us = User.objects.all()
us.delete()
```