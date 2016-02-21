define(['jquery'], function($) {


    var Interact = function(options) {
        this.config = $.extend({
            debug: false
        }, options);
    };

    Interact.prototype.log = function() {
        if (this.config.debug) {
            console.log('[Interact] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }

        return this;
    };

    Interact.prototype.init = function() {
        this.log('init()');
        var that = this,
            count = 0,
            gaInterval = setInterval(function() {
                that.log('check ga', count, typeof ga);
                if (typeof ga == "function") {
                    that.initHandlers();
                    clearInterval(gaInterval);
                    return;
                }
                count++;
                if (count > 20) { // Stop check after 10 seconds
                    clearInterval(gaInterval);
                }
            }, 500);

        return this;
    };

    Interact.prototype.initHandlers = function() {
        this.log('initHandlers()');

        var that = this,
            $body = $('body'),
            category = 'control';

        $body
            .on('etalia.control.search.change', function(e, data) {
                var obj = {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'search',
                    eventLabel:    data.expression
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.control.timespan.change', function(e, data) {
                var obj = {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'timespan',
                    eventLabel:    data.label,
                    eventValue:    data.value
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.control.cluster.change', function(e, data) {
                var obj = {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'cluster',
                    eventLabel:    data.label,
                    eventValue:    data.value
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.control.pinned.change', function(e, data) {
                var obj = {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   data.active ? 'pin' : 'unpin'
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.control.filters.change', function(e, data) {
                var obj = {
                    hitType: 'event',
                    eventCategory: category,
                    eventAction: (data.active ? 'select_' : 'deselect_') + data.group,
                    eventLabel: data.label,
                    eventValue: data.value
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });

        return this;
    };

    return new Interact().init();
});
