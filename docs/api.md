Publication API
===============

## Pin publication

url: ```[POST] /user/paper/pin```

**Request**

```json
{
    'pk':     (int),
    'source': (string),
}
```

**Response**

```json
{
    'is_liked':      (bool), # Pinned or not
    'is_ticked':     (bool), # Banned or not 
    'likes_counter': (int),  # ?
}
```

**Event**

```etalia.publication.pin```

```json
{
    'id':            (int),
    'is_liked':      (bool), # Pinned or not
    'is_ticked':     (bool), # Banned or not 
    'likes_counter': (int),  # ?
}
```

## Ban publication

url: ```[POST] /user/paper/ban```

**Request**

```json
{
    'pk':     (int),
    'source': (string),
}
```

**Response**

```json
{
    'is_liked':      (bool), # Pinned or not
    'is_ticked':     (bool), # Banned or not 
    'likes_counter': (int),  # ?
}
```

**Event**

```etalia.publication.ban```       

```json
{
    'id':            (int),
    'is_liked':      (bool), # Pinned or not
    'is_ticked':     (bool), # Banned or not 
    'likes_counter': (int),  # ?
}
```


## Add publication (to library)

url: ```[POST] /user/paper/add```

**Request**

```json
{
    'pk': (int),
}
```

**Response**

```json
{
    'success':         true,
    'trash_counter':   (int),
    'library_counter': (int),
    'likes_counter':   (int),
    'message':         (string)
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
    'trash_counter':   (int),
    'library_counter': (int),
    'likes_counter':   (int),
    'message':         (string)
}
```


## Trash publication (from library)

url: ```[POST] /user/paper/trash```

**Request**

```json
{
    'pk': (int),
}
```

**Response**

```json
{
    'success':         true,
    'trash_counter':   (int),
    'library_counter': (int),
    'likes_counter':   (int),
    'message':         (string)
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
    'trash_counter':   (int),
    'library_counter': (int),
    'likes_counter':   (int),
    'message':         (string)
}
```


## Restore publication (into library)

url: ```[POST] /user/paper/restore```

**Request**

```json
{
    'pk': (int),
}
```

**Response**

```json
{
    'success':         true,
    'trash_counter':   (int),
    'library_counter': (int),
    'likes_counter':   (int),
    'message':         (string)
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
    'trash_counter':   (int),
    'library_counter': (int),
    'likes_counter':   (int),
    'message':         (string)
}
```
