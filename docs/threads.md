THREADS
========

## Thread

Url: ```[GET] /threads/<pk>/```

**Response**

HTML Page


## Thread

Url: ```[POST] /threads/```

**Parameters**
```
{
    'type': enum                # required (see THREAD_TYPES in threads/constant.py)
    'privacy': enum             # required (see THREAD_PRIVACIES in threads/constant.py)
    'title': (str)              # required
    'paper': (int) | None
}
```

Url: ```[POST] /threads/<pk>/```

**Parameters**
```
{
    'title': (str)              # required
    'content': (str)
}
```

Url: ```[GET] /threads/<pk>/```
**Parameters**
```
```

**All Response**

```200 OK```
```
{
    "results": {
        "id": (int),
        "type": (enum),
        "title": (str),
        "content": (str),
        "owner": {
            "id": (int),
            "email": (string), 
            "first_name": (string),
            "last_name": (string),
            "url": (string),
            "photo_url": (string)
        },
        "privacy": (enum),
        "state": {
            "id": (int),
            "is_pinned": (bool),
            "is_banned": (bool),
            "is_joined": (bool),
            "is_left": (bool),
            "first_joined_at": (string),
            "last_left_at": (string),
        },
        "members": [
            {
                "id": (int),
                "email": (string), 
                "first_name": (string),
                "last_name": (string),
                "url": (string),
                "photo_url": (string)
            },
            ...
        ],
        "paper": {
            "id": (int), 
            "title": (str),
            "journal": {
                "id": (int),
                "title": (str),
                "short_title": (str),
                "url": (str)
            },
            "authors": [
                {
                    "id": (int), 
                    "first_name": (string),
                    "last_name": (string)
                },
            ...
            ], 
            "id_doi": (str),
            "id_pmi": (str),
            "id_arx": (str),
            "url": (str),
        },
        "posts": [
            {
                "id": (int),
                "thread": (int),
                "content": (string),
                "created": (string),
                "modified": (string),
                "author": {
                    "id": (int),
                    "email": (string), 
                    "first_name": (string),
                    "last_name": (string),
                    "url": (string),
                    "photo_url": (string)
                },
                "comments:" : [
                    {                                    
                        "id": (int),
                        "content": (string),
                        "created": (string),
                        "modified": (string),                
                        "author": {
                            "id": (int),
                            "email": (string), 
                            "first_name": (string),
                            "last_name": (string),
                            "url": (string),
                            "photo_url": (string)
                        },
                        "post": (int)
                    },
                    ...
                ]
            },
            ...
        ]
    }
}
```


## Post (new, update)

Url: ```[POST] /threads/<pk>/posts/```
Url: ```[POST] /threads/posts/<pk>/update```

**Parameters**
```
{
    'content':  (str)
}
```

**Response**

```200 OK```
```
{
    "results": {   
        "id": (int),
        "thread": (int),
        "content": (string),
        "created": (string),
        "modified": (string),
        "author": {
            "id": (int),
            "email": (string), 
            "first_name": (string),
            "last_name": (string),
            "url": (string),
        },
        "comments:" : [
            {                                    
                "id": (int),
                "content": (string),
                "created": (string),
                "modified": (string),                
                "author": {
                    "id": (int),
                    "email": (string), 
                    "first_name": (string),
                    "last_name": (string),
                    "url": (string),
                },
                "post": (int)
            },
            ...
        ]
    }
}
```

## Delete post

Url: ```[POST] /threads/posts/<pk>/delete```

**Response**
``` 200 OK```

## Comment (new, update)

Url: ```[POST] /threads/post/<pk>/comments/new```
Url: ```[POST] /threads/comments/<pk>/update```

**Parameters**
```
{
    'content':  (str)
}
```


**Response**

```200 OK```
```
{
    "results": {   
        "id": (int),
        "post": (int),
        "content": (string),
        "created": (string),
        "modified": (string),
        "author": {
            "id": (int),
            "email": (string), 
            "first_name": (string),
            "last_name": (string),
            "url": (string),
            "photo_url": (string)
        }
    }
}
```

## Delete comment

Url: ```[POST] /threads/comments/<pk>/delete```

**Response**
``` 200 OK```


## ThreadUserState

Url:  ```[POST] /threads/state```

**Request**
```json
{
    'action': (str),        # either 'pin', 'ban', 'join', 'leave'
}
```

**Response**
```json
{
    "results": {
        "id": (int),
        "is_pinned": (bool),
        "is_banned": (bool),
        "is_joined": (bool),
        "is_left": (bool)
    }
}
```