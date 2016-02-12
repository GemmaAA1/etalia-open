Publication API
===============

## Pin publication

url: ```[POST] /user/paper/pin```

**Request**

```json
{
    'id':     (int),
    'source': (string),
}
```

**Response**

```json
{
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```

**Event**

```etalia.publication.pin```

```json
{
    'id':            (int),
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':      (int),
        'ban':      (int),
        'trash':    (int),
        'library':  (int),
    }
}
```

## Ban publication

url: ```[POST] /user/paper/ban```

**Request**

```json
{
    'id':     (int),
    'source': (string),
}
```

**Response**

```json
{
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```

**Event**

```etalia.publication.ban```       

```json
{
    'id':            (int),
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```


## Add publication (to library)

url: ```[POST] /user/paper/add```

**Request**

```json
{
    'id': (int),
}
```

**Response**

```json
{
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
// OR
{
    'success': false,
    'message': (str),
}
```
           
**Event**

```etalia.publication.add```

```json
{
    'id':              (int),
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```


## Trash publication (from library)

url: ```[POST] /user/paper/trash```

**Request**

```json
{
    'id': (int),
}
```

**Response**

```json
{
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
// OR
{
    'success': false,
    'message': (str),
}
```
               
**Event**

```etalia.publication.trash```

```json
{
    'id':              (int),
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```


## Restore publication (into library)

url: ```[POST] /user/paper/restore```

**Request**

```json
{
    'id': (int),
}
```

**Response**

```json
{
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
// OR
{
    'success': false,
    'message': (str),
}
```
                      
**Event**

```etalia.publication.restore```

```json
{
    'id':              (int),
    'state': {
        'is_pinned':    (bool), # Pinned or not
        'is_banned':    (bool), # Banned or not
        'is_trashed':   (bool), # In user trash or not
        'is_added':     (bool), # In user library or not
    }
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```
