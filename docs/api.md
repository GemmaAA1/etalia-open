Publication API
===============

## TODO

* use http code 200 on success.
* use http code 500 on failure (with error message).
* remove 'success' and 'message' data.

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

```etalia.publication.pin``` (with Result object)


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

```etalia.publication.ban``` (with Result object)       


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

```etalia.publication.add``` (with Result object)


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

```etalia.publication.trash``` (with Result object)


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

```etalia.publication.restore``` (with Result object)     


## Empty user library trash

url: ```[POST] /user/library/trash/empty```

**Response**

```json
{
    'counter': {
        'pin':          (int),
        'ban':          (int),
        'trash':        (int),
        'library':      (int),
    }
}
```

**Event**

```etalia.publication.trash-clear```


## Event result object

```javascript
$('body').on('etalia.publication.pin', function(e, result) {

    // Paper id
    result.getId();           // (int)

    // Paper states
    result.isPinned();        // (bool)
    result.isBanned();        // (bool)
    result.isTrashed();       // (bool)
    result.isAdded();         // (bool)

    // Counters
    result.getPinCount();     // (int)
    result.getBanCount();     // (int)
    result.getTrashCount();   // (int)
    result.getLibraryCount(); // (int)
});

```
