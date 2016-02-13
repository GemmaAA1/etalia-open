## Paper neighbors publication

url: ```[XHR GET] /library/paper/neighbors```

**Request**

```json
{
    'id':       (int),
    'time-span: (int),      # in [1-365]
}
```

**Response**

```xml
<?xml encoding="UTF-8">
<response>
   <data>
        {
            'id':           (int)
            'time-span':    (int)
        }
   </data>
   <thumb-list>
    <![CDATA[
        HTML block
    ]]>
   </thumb-list>
</response>
```