define(['jquery', 'app/api'], function($, Api) {

    var List = function(options) {

        this.config = $.extend({
            debug: false,
            element: '#list',
            container: '.thumb-list'
        }, options);

        this.$element = $(this.config.element);

        this.loadXhr = null;
    };

    List.prototype.log = function() {
        if (this.config.debug) {
            console.log('[List] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }

        return this;
    };

    List.prototype.init = function() {
        this.log('init', this.$element);

        this.$element
            .on('click', '.thumb-pin', function(e) {
                var $thumb = $(e.target).closest('.thumb');

                Api.pin($thumb.data('id'), window.location.pathname);

                e.stopPropagation();
                e.preventDefault();
                return false;
            })
            .on('click', '.thumb-ban', function(e) {
                var $thumb = $(e.target).closest('.thumb');

                Api.ban($thumb.data('id'), window.location.pathname);

                e.stopPropagation();
                e.preventDefault();
                return false;
            })
            .on('click', '.thumb-library-add', function(e) {
                var $thumb = $(e.target).closest('.thumb');

                Api.add($thumb.data('id'));

                e.stopPropagation();
                e.preventDefault();
                return false;
            })
            .on('click', '.thumb-library-trash', function(e) {
                var $thumb = $(e.target).closest('.thumb');

                Api.trash($thumb.data('id'));

                e.stopPropagation();
                e.preventDefault();
                return false;
            })
            .on('click', '.thumb-library-restore', function(e) {
                var $thumb = $(e.target).closest('.thumb');

                Api.restore($thumb.data('id'));

                e.stopPropagation();
                e.preventDefault();
                return false;
            });

        return this;
    };

    List.prototype.load = function(controlsStates) {
        controlsStates = controlsStates || {};

        if (this.loadXhr) {
            this.loadXhr.abort();
        }

        var that = this,
            $body = $('body');

        this.loadXhr = $.ajax({
            method:   'GET',
            url:      that.$element.data('load-url'),
            data:     {'data': JSON.stringify(controlsStates)},
            dataType: 'xml'
        })
        .done(function(xml) {
            var eventData = {};

            var data = $(xml).find('data');
            if (data.length) {
                eventData.controlsStates = JSON.parse(data.text());

                // List title
                if (data.hasOwnProperty('number_of_papers')) {
                    that.$element.find('.list-title span').html(data['number_of_papers']);
                }
            }

            var list = $(xml).find('thumb-list');
            if (list.length) {
                that.$element.find(that.config.container).html($(list.text()));
            }

            $body.trigger('etalia.list.load', eventData);
        })
        .fail(function(xrh, status, error) {
            that.log('Load failure', xrh, status, error);
        });


        this.log('load', that.$element.data('load-url'));

        return this;
    };

    return List;
});

