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

**GA**

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'search', 
eventLabel:    expression,
eventValue:    none,
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

**GA**

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'timespan', 
eventLabel:    label,
eventValue:    value,
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

**GA**

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'cluster', 
eventLabel:    label,
eventValue:    value,
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

**GA**

active == True:

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'pin', 
eventLabel:    '',
eventValue:    none,
```

active == False:

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'unpin', 
eventLabel:    '',
eventValue:    none,
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

**GA**

if active == True

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'select_' + group,   // 'select_journal' | 'select_author'
eventLabel:    label,               // The filter title
eventValue:    value,               // The filter id
```
  
if active == False

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'deselect_' + group, // 'deselect_journal' | 'deselect_author'
eventLabel:    label,               // The filter title
eventValue:    value,               // The filter id
```  

## More Filters
     
**Event**

```etalia.control.more_filters.change```
         
**Data**

```json
{
    'group':  (string),      // The group label
}
```

**GA**

```
hitType:       'event',
eventCategory: 'control',
eventAction:   'more_filter,
eventLabel:    group,
eventValue:    none,           
```  

## Pin from thumb

**Event**

```etalia.control.thumb_pin.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id
    'active': (bool)
}
```

**GA**

if active == True:

```
hitType:       'event',
eventCategory: 'thumb',
eventAction:   'pin',
eventLabel:    label,
eventValue:    value,           
```

if active == False:

```
hitType:       'event',
eventCategory: 'thumb',
eventAction:   'unpin',
eventLabel:    label,
eventValue:    value,
```

## Ban from thumb

**Event**

```etalia.control.thumb_ban.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id
    'active': (bool)         // for consistency with pin (but cannot be set to False)
}
```

**GA**

```
hitType:       'event',
eventCategory: 'thumb',
eventAction:   'ban',
eventLabel:    label,
eventValue:    value,           
```

## Pin from detail

**Event**

```etalia.control.detail_pin.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id
    'active': (bool)
}
```

**GA**

if active == True:

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'pin',
eventLabel:    label,
eventValue:    value,           
```

if active == False:

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'unpin',
eventLabel:    label,
eventValue:    value,
```

## Ban from detail

**Event**

```etalia.control.detail_ban.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id
    'active': (bool)         // for consistency with pin (but cannot be set to False)
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'ban',
eventLabel:    label,
eventValue:    value,           
```

## Add from detail

**Event**

```etalia.control.detail_add.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id    
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'add',
eventLabel:    label,
eventValue:    value,  
```

## Trash from detail

**Event**

```etalia.control.detail_trash.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id    
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'trash',
eventLabel:    label,
eventValue:    value,  
```

## Restore from detail

**Event**

```etalia.control.detail_restore.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id     
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'restore',
eventLabel:    label,
eventValue:    value,  
```

## Tweet from detail

**Event**

```etalia.control.detail_tweet.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id     
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'tweet',
eventLabel:    label,
eventValue:    value,  
```

## Google+ from detail

**Event**

```etalia.control.detail_gplus.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id     
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'gplus',
eventLabel:    label,
eventValue:    value,  
```

## Email from detail

**Event**

```etalia.control.detail_email.change```
         
**Data**

```json
{
    'label':  (string),      // The paper title
    'value':  (int),         // The paper id     
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'email',
eventLabel:    label,
eventValue:    value,  
```

## Similar paper timespan from detail

**Event**

```etalia.control.detail_similar_timespan.change```
         
**Data**

```json
{
    'label':  (string),
    'value':  (int), 
}
```

**GA**

```
hitType:       'event',
eventCategory: 'detail',
eventAction:   'similar_timespan',
eventLabel:    label,
eventValue:    value,  
```



