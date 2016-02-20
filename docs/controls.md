Controls events
===============

## Search
     
**Event**

```etalia.control.search.change```
         
**Data**

```json
{
    'expression': (string)
}
```


## Timespan
     
**Event**

```etalia.control.timespan.change```
         
**Data**

```json
{
    'value': (int),
    'label': (string)
}
```


## Cluster
     
**Event**

```etalia.control.cluster.change```
         
**Data**

```json
{
    'value': (int),
    'label': (string)
}
```


## Pinned
     
**Event**

```etalia.control.pinned.change```
         
**Data**

```json
{
    'active': (bool)
}
```


## Filters
     
**Event**

```etalia.control.filters.change```
         
**Data**

```json
{
    'value':  (int),         // The filter id
    'label':  (string),      // The filter title
    'group':  (string),      // The group id/label
    'active': (bool)         // The filter state
}
```
