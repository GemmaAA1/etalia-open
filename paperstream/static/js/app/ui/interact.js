define(['jquery'], function($) {

    var $body = $('body'),
        category = 'Control';

    function initHandlers() {
        $body
            .on('etalia.control.search.change', function(e, data) {
                ga('send', {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'search',
                    eventLabel:    '',
                    eventValue:    ''
                });
            })
            .on('etalia.control.timespan.change', function(e, data) {
                ga('send', {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'timespan',
                    eventLabel:    data.label,
                    eventValue:    data.value
                });
            })
            .on('etalia.control.cluster.change', function(e, data) {
                ga('send', {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'cluster',
                    eventLabel:    data.label,
                    eventValue:    data.value
                });
            })
            .on('etalia.control.pinned.change', function(e, data) {
                ga('send', {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'pin',
                    eventLabel:    data.active ? 'on' : 'off',
                    eventValue:    data.active ? 1 : 0
                });
            })
            .on('etalia.control.filters.change', function(e, data) {
                if (data.active) {
                    ga('send', {
                        hitType: 'event',
                        eventCategory: category,
                        eventAction: 'filter.' + data.group,
                        eventLabel: data.label,
                        eventValue: data.value
                    });
                }
            });
    }

    // GA is not loaded by requireJs, we need to check if it's available
    var count = 0,
        gaInterval = setInterval(function() {
            if (typeof ga == "object") {
                initHandlers();
                clearInterval(gaInterval);
                return;
            }
            count++;
            if (count > 20) { // Stop check after 10 seconds
                clearInterval(gaInterval);
            }
        }, 500);
});
