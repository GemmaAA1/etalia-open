define(['app/app'], function(App) {

    var Tracking = function(options) {
        this.config = App.$.extend({
            debug: true
        }, options);
    };

    Tracking.prototype.log = function() {
        if (this.config.debug) {
            console.log('[Tracking] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }

        return this;
    };

    Tracking.prototype.init = function() {
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

    Tracking.prototype.initHandlers = function() {
        this.log('initHandlers()');

        var that = this;

        /**
         * Ui controls events
         */
        var category = 'control';
        App
            .on('etalia.control.search.change', function(data) {
                var obj = {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   'search',
                    eventLabel:    data.expression
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.control.timespan.change', function(data) {
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
            .on('etalia.control.cluster.change', function(data) {
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
            .on('etalia.control.pinned.change', function(data) {
                var obj = {
                    hitType:       'event',
                    eventCategory: category,
                    eventAction:   data.active ? 'pin' : 'unpin'
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })
            .on('etalia.control.filters.change', function(data) {
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
            .on('etalia.control.filters.more', function(data) {
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
         * Papers events
         */
        function sendPaperEvent(action, paper) {
            var obj = {
                hitType: 'event',
                eventCategory: 'paper', //result.getExtra('source'),
                eventAction: action,
                eventLabel: paper.get('title'),
                eventValue: paper.get('id')
            };
            that.log('ga.send()', obj);
            ga('send', obj);
        }
        App
            .on('etalia.paper.pin', function(paper) {
                sendPaperEvent('pin', paper);
            })
            .on('etalia.paper.unpin', function(paper) {
                sendPaperEvent('unpin', paper);
            })
            .on('etalia.paper.ban', function(paper) {
                sendPaperEvent('ban', paper);
            })
            .on('etalia.paper.unban', function(paper) {
                sendPaperEvent('unban', paper);
            })
            .on('etalia.paper.add', function(paper) {
                sendPaperEvent('add', paper);
            })
            .on('etalia.paper.trash', function(paper) {
                sendPaperEvent('trash', paper);
            })
            /*.on('etalia.paper.restore', function(paper) {
                sendPaperEvent('restore', paper);
            })*/
            /*.on('etalia.paper.trash-clear', function() {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'library',
                    eventAction: 'trash-clear'
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })*/
            .on('etalia.paper.share', function(paper, support) {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'paper', //data.source,
                    eventAction: 'share_' + support,
                    eventLabel: paper.get('title'),
                    eventValue: paper.get('id')
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });

        /**
         * Threads events
         */
        function sendThreadEvent(action, thread) {
            var obj = {
                hitType: 'event',
                eventCategory: 'thread', //result.getExtra('source'),
                eventAction: action,
                eventLabel: thread.get('title'),
                eventValue: thread.get('id')
            };
            that.log('ga.send()', obj);
            ga('send', obj);
        }
        App
            .on('etalia.thread.pin', function(thread) {
                sendThreadEvent('pin', thread);
            })
            .on('etalia.thread.unpin', function(thread) {
                sendThreadEvent('unpin', thread);
            })
            .on('etalia.thread.ban', function(thread) {
                sendThreadEvent('ban', thread);
            })
            .on('etalia.thread.unban', function(thread) {
                sendThreadEvent('unban', thread);
            })
            .on('etalia.thread.add', function(thread) {
                sendThreadEvent('add', thread);
            })
            .on('etalia.thread.trash', function(thread) {
                sendThreadEvent('trash', thread);
            })
            /*.on('etalia.thread.restore', function(thread) {
                sendThreadEvent('restore', thread);
            })*/
            /*.on('etalia.thread.trash-clear', function() {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'library',
                    eventAction: 'trash-clear'
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            })*/
            .on('etalia.thread.share', function(thread, support) {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'thread', //data.source,
                    eventAction: 'share_' + support,
                    eventLabel: thread.get('title'),
                    eventValue: thread.get('id')
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });

        /**
         * Detail events
         */
        /*App
            .on('etalia.detail.similar_timespan.change', function(data) {
                var obj = {
                    hitType: 'event',
                    eventCategory: 'detail',
                    eventAction: 'similar_timespan',
                    eventLabel: data.title,
                    eventValue: data.value
                };
                that.log('ga.send()', obj);
                ga('send', obj);
            });*/

        return this;
    };

    return new Tracking();
});
