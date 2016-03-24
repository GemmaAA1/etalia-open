THREADS
========

## Thread

Url: ```[GET] /threads/<pk>/```

**Response**

HTML Page


## Thread (new, update)

Url: ```[POST] /threads/new/```
Url: ```[POST] /threads/<pk>/update```

**Response**

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
        },
        "privacy": (enum),
        "members": [
            {
                "id": (int),
                "email": (string), 
                "first_name": (string),
                "last_name": (string),
                "url": (string),
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
                    "last_name": (string),
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
            },
            ...
        ],
    },
}
```


## Join thread

Url: ```[POST] /threads/<pk>/join```

**Response**

```200 OK```
```
{
}
```


## Post (new, update)

Url: ```[POST] /threads/<pk>/posts/new```
Url: ```[POST] /threads/posts/<pk>/update```

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
    },
}
```

## Delete post

Url: ```[POST] /threads/posts/<pk>/delete```

**Response**
``` 200 OK```

## Comment (new, update)

Url: ```[POST] /threads/post/<pk>/comments/new```
Url: ```[POST] /threads/comments/<pk>/update```

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
        },
    },
}
```

## Delete comment

Url: ```[POST] /threads/comments/<pk>/delete```

**Response**
``` 200 OK```
