define(['jquery'], function($) {


    var Interact = function(options) {
        this.config = $.extend({
            debug: true
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
            $body = $('body');

        /**
         * Ui controls events
         */
        var category = 'control';
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
            })
            .on('etalia.control.filters.more', function(e, data) {
                var obj = {
                    hitType: 'event',
                    eventCategory: category,
                    eventAction: 'more_filter',
                    eventLabel: data.group
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });

        /**
         * Publication events
         */
        $body
            .on('etalia.publication.pin', function(e, result) {
                var obj = {
                    hitType: 'event',
                    eventCategory: result.getExtra('source'),
                    eventAction: result.isPinned() ? 'pin' : 'unpin',
                    eventLabel: result.getExtra('title'),
                    eventValue: result.getId()
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.publication.ban', function(e, result) {
                var obj = {
                    hitType: 'event',
                    eventCategory: result.getExtra('source'),
                    eventAction: 'ban',
                    eventLabel: result.getExtra('title'),
                    eventValue: result.getId()
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.publication.add', function(e, result) {
                var obj = {
                    hitType: 'event',
                    eventCategory: result.getExtra('source'),
                    eventAction: 'add',
                    eventLabel: result.getExtra('title'),
                    eventValue: result.getId()
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.publication.trash', function(e, result) {
                var obj = {
                    hitType: 'event',
                    eventCategory: result.getExtra('source'),
                    eventAction: 'trash',
                    eventLabel: result.getExtra('title'),
                    eventValue: result.getId()
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.publication.restore', function(e, result) {
                var obj = {
                    hitType: 'event',
                    eventCategory: result.getExtra('source'),
                    eventAction: 'restore',
                    eventLabel: result.getExtra('title'),
                    eventValue: result.getId()
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.publication.trash-clear', function() {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'library',
                    eventAction: 'trash-clear'
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.publication.share', function(e, data) {
                var obj = {
                    hitType: 'event',
                    eventCategory: data.source,
                    eventAction: 'share_' + data.support,
                    eventLabel: data.title,
                    eventValue: data.id
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });

        /**
         * Detail events
         */
        $body.
            on('etalia.detail.similar_timespan.change', function(e, data) {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'detail',
                    eventAction: 'similar_timespan',
                    eventLabel: data.title,
                    eventValue: data.value
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });

        return this;
    };

    return new Interact().init();
});
