Threads API
==========



## Thread list

```[GET] /api/v1/threads```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| page          | (int)         | (optional) Number of the page (default = 1) | 
| pinned        | (int)         | (optional) Fetch only **pinned** threads if 1 (default = 0) |
| added         | (int)         | (optional) Fetch only **added** (to library) threads if 1 (default = 0) |
| trashed       | (int)         | (optional) Fetch only **trashed** (from library) threads if 1 (default = 0) |

Note: setting both added=1 AND trashed=1 should return a ```400 Bad request``` error.


**Response**

```200 OK```

```json
{
    "count":  (int),
    "next": (string|null),
    "prev": (string|null),
    "results": [
        {   
            "url": (string),                                    
            "id": (int),
            "type": (string)
            "title": (string),
            "content": (string),                    # @TODO truncate/size
            "comment_count": (int),
            "expert_count": (int),
            "created_at": (string),
            "updated_at": (string),
            "state": {
                "url": (string),
                "id": (int),
                "pinned": (bool),
                "banned": (bool),
                "added": (bool),
                "trashed": (bool)
            },
            "author": {
                "url": (string),
                "id": (int),
                "username": (string),               # @TODO username or first/last name
                "first_name": (string),
                "last_name": (string),
                "photo": (string)
            },
            "paper": {                              # @TODO could be not set
                "url": (string),
                "id": (int),
                "title": (string),
                "author": (string),                 # @TODO or author first/last name
            }
        },
        ...
    ]
}
```



## Thread create

```[POST] /api/v1/threads```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| type          | (int)         | (required) The thread type identifier | 
| author        | (int)         | (required) The thread author identifier | 
| title         | (string)      | (required) The thread title | 
| content       | (string)      | (required) The thread content | 



**Response**

```201 Created```

```json
{
    "url": (string),
    "id": (int),
    "type": (string)
    "title": (string),
    "content": (string),
    "comment_count": (int),
    "expert_count": (int),
    "created_at": (string),
    "updated_at": (string),
    "state": {
        "url": (string),
        "id": (int),
        "pinned": (bool),
        "banned": (bool),
        "added": (bool),
        "trashed": (bool)
    },
    "author": {
        "url": (string),
        "id": (int),
        "username": (string),               # @TODO username or first/last name
        "first_name": (string),
        "last_name": (string),
        "photo": (string)
    },
    "paper": {                              # @TODO could be not set
        "url": (string),
        "id": (int),
        "title": (string),
        "author": (string),                 # @TODO or author first/last name
    },
    "assets": [
        {
            "url": (string),                # @TODO ?
            "id": (int),
            "type": (string),
            "data": (string)                # @TODO Unexploitable ? Or parse regarding to type => data must be json.
        },
        ...
    ]
}
```

OR

```400 Bad Request```

@TODO review error format

```json
{
    "message": "Validation failed",    
    "errors": [
        {   
            "field": (string),                                    
            "error": (string)                                                
        },
        ...
    ]
}
```



## Thread detail

url: ```[GET] /api/v1/threads/{id}```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| id            | (int)         | (required) The thread identifier |


**Response**

```200 OK```

```json
{
    "url": (string),
    "id": (int),
    "type": (string)
    "title": (string),
    "content": (string),
    "comment_count": (int),
    "expert_count": (int),
    "created_at": (string),
    "updated_at": (string),
    "state": {
        "url": (string),
        "id": (int),
        "pinned": (bool),
        "banned": (bool),
        "added": (bool),
        "trashed": (bool)
    },
    "author": {
        "url": (string),
        "id": (int),
        "username": (string),               # @TODO username or first/last name
        "first_name": (string),
        "last_name": (string),
        "photo": (string)
    },
    "paper": {                              # @TODO could be not set
        "url": (string),
        "id": (int),
        "title": (string),
        "author": (string),                 # @TODO or author first/last name
    },
    "assets": [
        {
            "url": (string),                # @TODO ?
            "id": (int),
            "type": (string),
            "data": (string)                # @TODO Unexploitable ? Or parse regarding to type => data must be json.
        },
        ...
    ]
}
```
                   
OR

```404 Not found```



## Thread update

```[PUT] /api/v1/threads/{id}```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| id            | (int)         | (required) The thread identifier | 
| type          | (int)         | (required) The thread type identifier | 
| author        | (int)         | (required) The thread author identifier | 
| title         | (string)      | (required) The thread title | 
| content       | (string)      | (required) The thread content | 



**Response**

```200 OK```

```json
{
    "url": (string),
    "id": (int),
    "type": (string)
    "title": (string),
    "content": (string),
    "comment_count": (int),
    "expert_count": (int),
    "created_at": (string),
    "updated_at": (string),
    "state": {
        "url": (string),
        "id": (int),
        "pinned": (bool),
        "banned": (bool),
        "added": (bool),
        "trashed": (bool)
    },
    "author": {
        "url": (string),
        "id": (int),
        "username": (string),               # @TODO username or first/last name
        "first_name": (string),
        "last_name": (string),
        "photo": (string)
    },
    "paper": {                              # @TODO could be not set
        "url": (string),
        "id": (int),
        "title": (string),
        "author": (string),                 # @TODO or author first/last name
    },
    "assets": [
        {
            "url": (string),                # @TODO ?
            "id": (int),
            "type": (string),
            "data": (string)                # @TODO Unexploitable ? Or parse regarding to type => data must be json.
        },
        ...
    ]
}
```

OR

```400 Bad Request```

@TODO review error format

```json
{
    "message": "Validation failed",    
    "errors": [
        {   
            "field": (string),                                    
            "error": (string)                                                
        },
        ...
    ]
}
```



## Thread partial update

```[PATCH] /api/v1/threads/{id}```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| id            | (int)         | (required) The thread identifier | 
| type          | (int)         | (optional) The thread type identifier | 
| author        | (int)         | (optional) The thread author identifier | 
| title         | (string)      | (optional) The thread title | 
| content       | (string)      | (optional) The thread content | 



**Response**

```200 OK```

```json
{
    "url": (string),
    "id": (int),
    "type": (string)
    "title": (string),
    "content": (string),
    "comment_count": (int),
    "expert_count": (int),
    "created_at": (string),
    "updated_at": (string),
    "state": {
        "url": (string),
        "id": (int),
        "pinned": (bool),
        "banned": (bool),
        "added": (bool),
        "trashed": (bool)
    },
    "author": {
        "url": (string),
        "id": (int),
        "username": (string),               # @TODO username or first/last name
        "first_name": (string),
        "last_name": (string),
        "photo": (string)
    },
    "paper": {                              # @TODO could be not set
        "url": (string),
        "id": (int),
        "title": (string),
        "author": (string),                 # @TODO or author first/last name
    },
    "assets": [
        {
            "url": (string),                # @TODO ?
            "id": (int),
            "type": (string),
            "data": (string)                # @TODO Unexploitable ? Or parse regarding to type => data must be json.
        },
        ...
    ]
}
```

OR

```400 Bad Request```

@TODO review error format

```json
{
    "message": "Validation failed",    
    "errors": [
        {   
            "field": (string),                                    
            "error": (string)                                                
        },
        ...
    ]
}
```




## Thread neighbors list

url: ```[GET] /api/v1/threads/{id}/neighbors```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| id            | (int)         | (required) The thread identifier |


**Response**

```200 OK```

```json
{
    "count":  (int),
    "next": (string|null),
    "prev": (string|null),
    "results": [
        {   
            "url": (string),                                    
            "id": (int),
            "type": (string)
            "title": (string),
            "content": (string),                    # @TODO truncate/size
            "comment_count": (int),
            "expert_count": (int),
            "created_at": (string),
            "updated_at": (string),
            "state": {
                "url": (string),
                "id": (int),
                "pinned": (bool),
                "banned": (bool),
                "added": (bool),
                "trashed": (bool)
            },
            "author": {
                "url": (string),
                "id": (int),
                "username": (string),               # @TODO username or first/last name
                "first_name": (string),
                "last_name": (string),
                "photo": (string)
            },
            "paper": {
                "url": (string),
                "id": (int),
                "title": (string),
                "author": (string),                 # @TODO or author first/last name
            }
        },
        ...
    ]
}
```



## Thread answers list

url: ```[GET] /api/v1/threads/{id}/answers```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| id            | (int)         | (required) The thread identifier |


**Response**

```200 OK```

```json
{
    "count":  (int),
    "next": (string|null),
    "prev": (string|null),
    "results": [
        {   
            "url": (string),                                    
            "id": (int),
            "thread": (string),                      # FQ Url
            "content": (string),
            "created_at": (string),
            "updated_at": (string),
            "author": {
                "url": (string),
                "id": (int),
                "username": (string),               # @TODO username or first/last name
                "first_name": (string),
                "last_name": (string),
                "photo": (string)
            },
            "comments:" : [
                {
                    "url": (string),                                    
                    "id": (int),
                    "content": (string),
                    "author": {
                        "url": (string),
                        "id": (int),
                        "username": (string),               # @TODO username or first/last name
                        "first_name": (string),
                        "last_name": (string)
                    },
                },
                ...
            ]
        },
        ...
    ]
}
```



## Thread answers create

url: ```[POST] /api/v1/threads/{id}/answers```

**Parameters**

| Name          | Type          | Description   |
| ------------- | ------------- | ------------- |
| id            | (int)         | (required) The thread identifier |


**Response**

```200 OK```

```json
{
    "count":  (int),
    "next": (string|null),
    "prev": (string|null),
    "results": [
        {   
            "url": (string),                                    
            "id": (int),
            "thread": (string),                      # FQ Url
            "content": (string),
            "created_at": (string),
            "updated_at": (string),
            "author": {
                "url": (string),
                "id": (int),
                "username": (string),               # @TODO username or first/last name
                "first_name": (string),
                "last_name": (string),
                "photo": (string)
            },
            "comments:" : [
                {
                    "url": (string),                                    
                    "id": (int),
                    "content": (string),
                    "author": {
                        "url": (string),
                        "id": (int),
                        "username": (string),               # @TODO username or first/last name
                        "first_name": (string),
                        "last_name": (string)
                    },
                },
                ...
            ]
        },
        ...
    ]
}
```

