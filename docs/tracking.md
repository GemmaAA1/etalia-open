Interaction tracking
====================

Google analytics custom events.

                                                                  
## Search control
     
On ```etalia.control.search.change```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'control',
    eventAction:   'search', 
    eventLabel:    (string), // The search expression
});
```


## Timespan control
     
On ```etalia.control.timespan.change```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'control',
    eventAction:   'timespan', 
    eventLabel:    (string),    // The selected timespan label
    eventValue:    (int),       // The selected timespan value
});
```


## Cluster control
     
On ```etalia.control.cluster.change```
   
```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'control',
    eventAction:   'cluster', 
    eventLabel:    (string),    // The selected cluster color
    eventValue:    (int),       // The selected cluster value
});
```


## Pinned control
     
On ```etalia.control.pinned.change```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'control',
    eventAction:   (string),  // 'pin' or 'unpin'
}); 
```

TODO: action = 'pin', value = 1 | 0 ?


## Filters control
     
On ```etalia.control.filters.change```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'control',
    eventAction:   (string),    // if 'active 'select_' + group_name, else 'deselect_' + group_name   (ex: 'select_journal')
    eventLabel:    (string),    // The filter title
    eventValue:    (int),       // The filter id
}); 
```
     
On ```etalia.control.filters.more```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'control',
    eventAction:   'more_filter,
    eventLabel:    (string),        // The group name
});            
```  


## Publication

On ```etalia.publication.pin```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: (string),    // The event emiter : 'thumb' or 'detail'
    eventAction:   (string),    // The publication pin state ('pin' or 'unpin')
    eventLabel:    (string),    // The publication title
    eventValue:    (int),       // The publication id
});                       
```

On ```etalia.publication.ban```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: (string),    // The event emiter : 'thumb' or 'detail'
    eventAction:   'ban',    
    eventLabel:    (string),    // The publication title
    eventValue:    (int),       // The publication id
});                       
```

On ```etalia.publication.add```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: (string),    // The event emiter : 'thumb' or 'detail'
    eventAction:   'add',    
    eventLabel:    (string),    // The publication title
    eventValue:    (int),       // The publication id
});                       
```

On ```etalia.publication.trash```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: (string),    // The event emiter : 'thumb' or 'detail'
    eventAction:   'trash',    
    eventLabel:    (string),    // The publication title
    eventValue:    (int),       // The publication id
});                       
```

On ```etalia.publication.restore```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: (string),    // The event emiter : 'thumb' or 'detail'
    eventAction:   'restore',    
    eventLabel:    (string),    // The publication title
    eventValue:    (int),       // The publication id
});                       
```

On ```etalia.publication.share```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: (string),    // The event emiter : 'thumb' or 'detail'
    eventAction:   (string),    // 'share_' + support ('email', 'google-plus' or 'twitter')    
    eventLabel:    (string),    // The publication title
    eventValue:    (int),       // The publication id
});                       
```

On ```etalia.publication.trash-clear```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'library,    
    eventAction:   'trash-clear'
});                       
```

## Paper detail

On ```etalia.detail.similar_timespan.change```

```javascript
ga('send', {
    hitType:       'event',
    eventCategory: 'detail',
    eventAction:   'similar_timespan',    
    eventLabel:    (string),    // The selected timespan label
    eventValue:    (string),    // The selected timespan value
}); 
```
