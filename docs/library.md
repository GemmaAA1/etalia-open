## Paper neighbors publication

**TODO** 

* rewrite to ```[XHR GET] /library/paper/{id}/neighbors``` ('rest' design)

url: ```[XHR GET] /library/paper/neighbors```

**Request**

**TODO** 

* remove 'id'
* rename 'time-span' to 'time_span' (to match feeds naming convention)

```json
{
    'id':       (int),
    'time-span: (int),      # in [1-365]
}
```

**Response**

**TODO** 

* remove 'id'
* rename 'time-span' to 'time_span' (to match feeds naming convention)

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
