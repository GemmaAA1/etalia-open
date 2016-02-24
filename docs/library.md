## Paper neighbors publication

**TODO** 

url: ```[XHR GET] /library/paper/{id}/neighbors```

**Request**

```json
{
    'time_span: (int),      # in [1-365]
}
```

**Response**

```xml
<?xml encoding="UTF-8">
<response>
   <data>
        {
            'time_span':    (int)
        }
   </data>
   <thumb-list>
    <![CDATA[
        HTML block
    ]]>
   </thumb-list>
</response>
```
