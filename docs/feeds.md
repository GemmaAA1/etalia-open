Publication Feeds
=================

<a name="stream"></a>
## Stream

Initial list of publications.

Url: ```[GET] /feed/stream/```

**Response**

    HTML page


<a name="stream-endless-scroll"></a>
#### Stream endless scroll 

(Ajax) Load next publications.

Url: ```[XHR GET] /feed/stream/```

**Request**

```json
{
    'journals':     (int)[],
    'authors':      (int)[],
    'cluster':      (int) or null,
    'pin':          (bool),
    'search_query': (str),
    // Paging ?
}
```

**Response**

    HTML block


<a name="stream-filtering"></a>
#### Stream filtering

(Ajax) Reset the list of publications.

Url: ```[XHR GET] /feed/stream/xml```

**Request**

```json
{
    'filters': [
        {
            'id':   (string),       # Values: journal, author
            'pk':   (int)[]
        },
        ...
    ],
    'time_span':    (int),          # Values: 7, 30, 60
    'cluster':      (int) or null,  # Values: null, 0, 1, 2, 3
    'pin':          (bool),
    'search_query': (str),
}
```

**Response**

```xml
<?xml encoding="UTF-8">
<response>
   <data>
        {
            'filters': [
                {
                    'id': (string),                    # Values: journal, author
                    'entries': [
                        {
                            'pk':         (int),
                            'name':       (string),
                            'count':      (int|null),
                            'is_checked': (bool)
                        },
                        ...
                    ]
                },
                ...
            ],
            'time_span':        (int),
            'cluster':          (int) or null,
            'pin':              (bool),
            'search_query':     (str),
            'number_of_papers': (int),
        }
   </data>
   <thumb-list>
    <![CDATA[
        HTML block
    ]]>
   </thumb-list>
</response>
```



<a name="trend"></a>
## Trend

Initial list of publications.

Url: ```[GET] /feed/trend/```

**Request** / **Response** same as [Stream](#stream)


<a name="trend-endless-scroll"></a>
#### Trend endless scroll 

Url: ```[XHR GET] /feed/trend/```

**Request** / **Response** same as [Stream endless scroll](#stream-endless-scroll)

                                 
#### Trend filtering
 
url: ```[XHR GET] /feed/trend/xml```

**Request** / **Response** same as [Stream filtering](#stream-filtering)
